#!/usr/bin/env python3
"""Volumetric parcellation construction for Chimera.

This module decomposes the original ~1,780-line ``Chimera.build_parcellation``
method into a :class:`ParcellationBuilder` whose named methods each handle one
step of the pipeline.

The original method relied heavily on ``"<name>" in locals()`` checks to carry
optional state between blocks.  That pattern cannot survive an "extract method"
refactor, so every such piece of cross-step state is now an explicit attribute
on the builder instance (initialised to ``None`` in :meth:`__init__`).  The
per-supra-region temporaries (``lh_supra_parc`` etc.) are reset at the start of
every iteration, which is equivalent to the ``del`` statements the original
executed at the end of each iteration.

The heavy per-block logic is transcribed verbatim from the original; each method
rebinds the read-only locals it needs (``deriv_dir``, ``fullid`` …) at its top so
the bodies remain a faithful copy.
"""

import copy
import os
import shutil
import subprocess
import uuid
from datetime import datetime
from glob import glob
from pathlib import Path

import clabtoolkit.bidstools as cltbids
import clabtoolkit.freesurfertools as cltfree
import clabtoolkit.imagetools as cltimg
import clabtoolkit.misctools as cltmisc
import clabtoolkit.parcellationtools as cltparc
import clabtoolkit.segmentationtools as cltseg
import nibabel as nib
import numpy as np
from clabtoolkit.colorstools import ColorTableLoader

from .config_manager import _pipeline_info
from .parcellation import _mix_side_prop, create_extra_regions_parc
from .processing import launch_fsl_first
from .region_table import build_lut_header


class ParcellationBuilder:
    """Build the Chimera parcellation for a single T1w image.

    Parameters mirror the original ``Chimera.build_parcellation`` signature.

    Parameters
    ----------
    chim : Chimera
        The configured Chimera object (already carrying ``parc_dict``,
        ``supra_dict`` and ``parc_code``, and with templates prepared).

    t1 : str
        T1-weighted image in the BIDs format.

    bids_dir : str
        BIDs dataset directory.

    deriv_dir : str
        BIDs derivative directory.

    fssubj_dir : str
        FreeSurfer subjects directory.

    growwm : str or int or list
        Grow of GM labels inside the white matter in mm.

    bool_mixwm : bool
        If True, the grown voxels mix with the cortical labels.

    force : bool
        Overwrite the results.

    pipe_dict : dict
        Pipeline configuration dictionary.  Resolved with :func:`_pipeline_info`
        when ``None``.
    """

    def __init__(
        self,
        chim,
        t1: str,
        bids_dir: str,
        deriv_dir: str = None,
        fssubj_dir: str = None,
        growwm: str | int | list = None,
        bool_mixwm: bool = False,
        force: bool = False,
        pipe_dict: dict = None,
    ):
        self.chim = chim
        self.t1 = t1
        self.bids_dir = bids_dir
        self.deriv_dir = deriv_dir
        self.fssubj_dir = fssubj_dir
        self.growwm = growwm
        self.bool_mixwm = bool_mixwm
        self.force = force
        self.pipe_dict = pipe_dict if pipe_dict is not None else _pipeline_info()

        self.chim_code = chim.parc_code

        # --- Path / identity state (set by _resolve_paths) -----------------
        self.anat_dir = None
        self.t1_name = None
        self.ent_dict = None
        self.ent_dict_fullid = None
        self.fullid = None
        self.path_cad = None
        self.out_dir = None
        self.prev_chims = None
        self.supra_names = None
        self.bool_ctx = False

        # --- FreeSurfer / base-segmentation state --------------------------
        self.lut_dict = None
        self.sub2proc = None
        self.cont_tech_freesurfer = None
        self.cont_image_freesurfer = None
        self.freesurfer_license_file = None
        self.cont_tech_ants = None
        self.cont_image_ants = None
        self.cont_tech_fsl = None
        self.cont_image_fsl = None
        self.aseg_parc = None
        self.extra_parc = None
        self.glob_header_info = None
        self.gm_sub_names = None
        self.t1_image = None
        self.affine = None
        self.files2del = None

        # --- Accumulating parcellations ------------------------------------
        self.lh_parc = None
        self.rh_parc = None
        self.mid_parc = None
        self.lh2refill = None
        self.rh2refill = None
        self.first_parc = None

        # --- Per-supra-region temporaries (reset each iteration) -----------
        self.lh_supra_parc = None
        self.rh_supra_parc = None
        self.mid_supra_parc = None

        # --- Region counts -------------------------------------------------
        self.nlh_subc = None
        self.nrh_subc = None
        self.nmid_subc = None

    # -- Convenience accessors so transcribed bodies read like the original --
    @property
    def parc_dict(self):
        return self.chim.parc_dict

    @property
    def supra_dict(self):
        return self.chim.supra_dict

    # ------------------------------------------------------------------ run
    def run(self):
        """Build the parcellation, computing it only if the outputs are missing."""
        self._resolve_paths()

        if not self._outputs_already_exist():
            self._load_freesurfer_lut()
            self._run_freesurfer_base()
            self._process_gm_supraregions()
            self._count_subcortical_regions()
            self._assemble_and_write()
            self._cleanup()

    # -------------------------------------------------------- step: paths
    def _resolve_paths(self):
        """Resolve BIDs entities and the Chimera output directory."""
        t1 = self.t1
        bids_dir = self.bids_dir
        deriv_dir = self.deriv_dir

        if not os.path.isfile(t1):
            raise ValueError("Please provide a valid T1 image.")

        # Getting the entities from the name
        self.anat_dir = os.path.dirname(t1)
        t1_name = os.path.basename(t1)
        self.t1_name = t1_name
        self.ent_dict = cltbids.str2entity(t1_name)

        temp_entities = t1_name.split("_")[:-1]
        fullid = "_".join(temp_entities)
        self.fullid = fullid
        self.ent_dict_fullid = cltbids.str2entity(fullid)

        if "ses" in self.ent_dict.keys():
            path_cad = (
                "sub-" + self.ent_dict["sub"] + os.path.sep + "ses-" + self.ent_dict["ses"]
            )
        else:
            path_cad = "sub-" + self.ent_dict["sub"]
        self.path_cad = path_cad

        # Creating Chimera directories
        if deriv_dir is None:
            out_dir = os.path.join(bids_dir, "derivatives", "chimera", path_cad, "anat")
        else:
            out_dir = os.path.join(deriv_dir, "chimera", path_cad, "anat")

        # Create the Chimera directory if it does not exist
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        self.out_dir = out_dir

        # Detecting all the previous parcellations generated for the subject
        self.prev_chims = glob(os.path.join(out_dir, "*chimera*_dseg.nii.gz"))

        supra_names = list(self.parc_dict.keys())
        # Detecting if Cortical is on the list of supra-regions
        self.bool_ctx = False
        if "Cortical" in supra_names:
            # Remove it from the list
            supra_names.remove("Cortical")
            self.bool_ctx = True
        self.supra_names = supra_names

    # --------------------------------------------- step: existence check
    def _outputs_already_exist(self) -> bool:
        """Return ``True`` when every requested output file is already present."""
        force = self.force
        chim_code = self.chim_code
        fullid = self.fullid
        chim_dir = self.out_dir

        if force:
            return False

        bool_chim_exist = True

        if self.bool_ctx:

            # Atributes for the cortical parcellation
            atlas_names = self.parc_dict["Cortical"]["parcels"]
            for at_name in atlas_names:
                ## -------- Cortical parcellation for the left hemisphere -------
                # Creating the name for the output file

                if self.growwm is None:
                    self.growwm = ["0"]

                for ngrow in np.arange(len(self.growwm)):
                    if self.growwm[ngrow] == "0":
                        out_vol_name = fullid + "_" + at_name + "_dseg.nii.gz"
                    else:

                        ent_dict = cltbids.str2entity(at_name)
                        if "desc" in ent_dict.keys():
                            if self.bool_mixwm:
                                ent_dict["desc"] = (
                                    ent_dict["desc"]
                                    + "grow"
                                    + str(self.growwm[ngrow])
                                    + "mm+mixwm"
                                )
                            else:
                                ent_dict["desc"] = (
                                    ent_dict["desc"]
                                    + "grow"
                                    + str(self.growwm[ngrow])
                                    + "mm"
                                )
                            tmp_str = cltbids.entity2str(ent_dict)
                            out_vol_name = fullid + "_" + tmp_str + "_dseg.nii.gz"
                        else:
                            if self.bool_mixwm:
                                out_vol_name = (
                                    fullid
                                    + "_"
                                    + at_name
                                    + "_desc-grow"
                                    + str(self.growwm[ngrow])
                                    + "mm+mixwm_dseg.nii.gz"
                                )
                            else:
                                out_vol_name = (
                                    fullid
                                    + "_"
                                    + at_name
                                    + "_desc-grow"
                                    + str(self.growwm[ngrow])
                                    + "mm_dseg.nii.gz"
                                )

                    chim_parc_name = cltbids.replace_entity_value(
                        out_vol_name, {"atlas": "chimera" + chim_code}
                    )
                    chim_parc_file = os.path.join(str(chim_dir), chim_parc_name)
                    chim_parc_lut = os.path.join(
                        str(chim_dir),
                        cltbids.replace_entity_value(
                            chim_parc_name, {"extension": "lut"}
                        ),
                    )
                    chim_parc_tsv = os.path.join(
                        str(chim_dir),
                        cltbids.replace_entity_value(
                            chim_parc_name, {"extension": "tsv"}
                        ),
                    )

                    if (
                        not os.path.isfile(chim_parc_file)
                        or not os.path.isfile(chim_parc_lut)
                        or not os.path.isfile(chim_parc_tsv)
                        or force
                    ):
                        bool_chim_exist = bool_chim_exist & False
                    else:
                        bool_chim_exist = bool_chim_exist & True

        else:
            out_vol_name = fullid + "_dseg.nii.gz"

            chim_parc_name = cltbids.insert_entity(
                out_vol_name, {"atlas": "chimera" + chim_code}
            )
            chim_parc_file = os.path.join(str(chim_dir), chim_parc_name)
            chim_parc_lut = os.path.join(
                str(chim_dir),
                cltbids.replace_entity_value(chim_parc_name, {"extension": "lut"}),
            )
            chim_parc_tsv = os.path.join(
                str(chim_dir),
                cltbids.replace_entity_value(chim_parc_name, {"extension": "tsv"}),
            )

            if (
                not os.path.isfile(chim_parc_file)
                or not os.path.isfile(chim_parc_lut)
                or not os.path.isfile(chim_parc_tsv)
                or force
            ):
                bool_chim_exist = bool_chim_exist & False
            else:
                bool_chim_exist = bool_chim_exist & True

        return bool(bool_chim_exist)

    # ------------------------------------------- step: FreeSurfer color LUT
    def _load_freesurfer_lut(self):
        """Read the FreeSurfer color LUT (locally or through a container)."""
        pipe_dict = self.pipe_dict

        if pipe_dict["packages"]["freesurfer"]["cont_tech"] != "local":
            cont_tech = pipe_dict["packages"]["freesurfer"]["cont_tech"]
            cont_image = pipe_dict["packages"]["freesurfer"]["container"]

            # Detecting the FreeSurfer home directory using the container
            cmd_bashargs = ["echo", "$FREESURFER_HOME"]
            cmd_cont = cltmisc.generate_container_command(
                cmd_bashargs, cont_tech, cont_image
            )  # Generating container command
            out_cmd = subprocess.run(
                cmd_cont, stdout=subprocess.PIPE, text=True
            )
            fslut_file_cont = os.path.join(
                out_cmd.stdout.split("\n")[0], "FreeSurferColorLUT.txt"
            )

            tmp_name = str(uuid.uuid4())
            cmd_bashargs = ["cp", "replace_cad", "/tmp/" + tmp_name]
            cmd_cont = cltmisc.generate_container_command(
                cmd_bashargs, cont_tech, cont_image
            )

            # Replace the element of the list equal to replace_cad by the path of the lut file
            cmd_cont = [w.replace("replace_cad", fslut_file_cont) for w in cmd_cont]
            subprocess.run(cmd_cont, stdout=subprocess.PIPE, text=True)
            fslut_file = os.path.join("/tmp", tmp_name)

            ######## ------------- Reading FreeSurfer color lut table ------------ #
            self.lut_dict = ColorTableLoader.read_luttable(fslut_file)
            os.remove(fslut_file)

        else:

            fshome_dir = os.getenv("FREESURFER_HOME")
            fslut_file = os.path.join(fshome_dir, "FreeSurferColorLUT.txt")

            ######## ------------- Reading FreeSurfer color lut table ------------ #
            self.lut_dict = ColorTableLoader.read_luttable(fslut_file)

    # ------------------------------------- step: FreeSurfer + base segments
    def _run_freesurfer_base(self):
        """Run FreeSurfer, build the aseg + extra parcellations and the buffers."""
        pipe_dict = self.pipe_dict
        t1 = self.t1
        force = self.force
        fssubj_dir = self.fssubj_dir
        fullid = self.fullid

        ######## ----- Running FreeSurfer if it was not previously computed ------ #
        sub2proc = cltfree.FreeSurferSubject(fullid, subjs_dir=fssubj_dir)
        self.sub2proc = sub2proc

        self.cont_tech_freesurfer = pipe_dict["packages"]["freesurfer"]["cont_tech"]
        self.cont_image_freesurfer = pipe_dict["packages"]["freesurfer"]["container"]
        self.freesurfer_license_file = pipe_dict["packages"]["freesurfer"]["license"]
        self.cont_tech_ants = pipe_dict["packages"]["ants"]["cont_tech"]
        self.cont_image_ants = pipe_dict["packages"]["ants"]["container"]
        self.cont_tech_fsl = pipe_dict["packages"]["fsl"]["cont_tech"]
        self.cont_image_fsl = pipe_dict["packages"]["fsl"]["container"]

        cont_tech_freesurfer = self.cont_tech_freesurfer
        cont_image_freesurfer = self.cont_image_freesurfer
        freesurfer_license_file = self.freesurfer_license_file

        # Running FreeSurfer if it was not previously computed is mandatory
        sub2proc.launch_freesurfer(
            force=force,
            t1w_img=t1,
            cont_tech=cont_tech_freesurfer,
            cont_image=cont_image_freesurfer,
            fs_license=freesurfer_license_file,
        )

        nii_image = os.path.join(
            sub2proc.subjs_dir, sub2proc.subj_id, "tmp", "aparc+aseg.nii.gz"
        )
        mgz_image = os.path.join(
            sub2proc.subjs_dir, sub2proc.subj_id, "mri", "aparc+aseg.mgz"
        )

        if not os.path.isfile(nii_image):
            sub2proc.conform2native(
                mgz_conform=mgz_image,
                nii_native=nii_image,
                cont_image=cont_image_freesurfer,
                cont_tech=cont_tech_freesurfer,
                force=force,
            )

        if self.aseg_parc is None:
            self.aseg_parc = cltparc.Parcellation(
                parc_file=nii_image, color_table=self.lut_dict
            )

            # Adjusting the index, color and names to contain only the FreeSurfer
            #  labels present values in the aparc+aseg parcellation
            self.aseg_parc.adjust_values()

        # Creating the parcellation for the extra regions
        self.extra_parc = create_extra_regions_parc(aparc=nii_image)

        # Remove the nifti file
        os.remove(nii_image)

        # Building the main header information for the LUT file
        self.glob_header_info = build_lut_header(self.chim)

        # Processing that will be perfomed for multiple supra-regions
        gm_sub_names = list(self.parc_dict.keys())

        # Remove 'Cortical', 'GyralWM' and 'WhiteMatter' from the list
        if "Cortical" in gm_sub_names:
            gm_sub_names.remove("Cortical")

        if "GyralWM" in gm_sub_names:
            gm_sub_names.remove("GyralWM")

        if "WhiteMatter" in self.supra_names:
            gm_sub_names.remove("WhiteMatter")
        self.gm_sub_names = gm_sub_names

        # Taking the image dimensions and the affine matrix in native space
        t1_image = nib.load(t1)
        self.t1_image = t1_image
        self.affine = t1_image.affine

        # Create a numpy array with the same dimensions of the T1 image and fill it
        # with zeros. The array will be used to store the parcellation.
        ref_image = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
        self.lh2refill = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
        self.rh2refill = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)

        # Creating the parcellation objects
        self.lh_parc = cltparc.Parcellation(
            parc_file=ref_image, affine=self.affine
        )  # It will include the parcellation for the left hemisphere
        self.rh_parc = copy.deepcopy(
            self.lh_parc
        )  # It will include the parcellation for the right hemisphere
        self.mid_parc = copy.deepcopy(
            self.lh_parc
        )  # It will include the parcellation for structures that do not belong to any hemisphere

        self.files2del = []  # Temporal files that will be deleted

    # ------------------------------------- step: per-supra-region processing
    def _process_gm_supraregions(self):
        """Process every gray-matter supra-region and accumulate the results."""
        for supra in self.gm_sub_names:
            self._process_one_supraregion(supra)

    def _process_one_supraregion(self, supra):
        # Reset the per-iteration temporaries (equivalent to the original's
        # ``del`` at the end of each iteration).
        self.lh_supra_parc = None
        self.rh_supra_parc = None
        self.mid_supra_parc = None

        # Getting the information of the common atributes
        atlas_code = self.parc_dict[supra]["code"]
        atlas_str = self.parc_dict[supra]["atlas"]
        atlas_desc = self.parc_dict[supra]["description"]
        atlas_cita = self.parc_dict[supra]["citation"]
        atlas_src = self.parc_dict[supra]["source"]
        atlas_ref = self.parc_dict[supra]["reference"]
        atlas_parcs = self.parc_dict[supra]["parcels"]
        atlas_mask = self.parc_dict[supra]["mask"]
        atlas_type = self.parc_dict[supra]["type"]

        # Check if there is a previous parcellation for the supra-region
        prev_parc = self._find_previous_parc(supra)

        if prev_parc is not None:
            self._load_supra_from_previous(supra, atlas_code, prev_parc)
        else:
            if atlas_type == "spam":
                # Verifying the existence of the threshold
                if "probthresh" not in self.parc_dict[supra].keys():
                    spam_thresh = 0.05
                else:
                    spam_thresh = self.parc_dict[supra]["probthresh"]
            else:
                spam_thresh = None

            deriv_fold = self.parc_dict[supra]["deriv_volfold"]
            proc_dict = self.parc_dict[supra]["processing"]

            if proc_dict["method"] == "comform2native":
                self._process_comform2native(
                    supra, atlas_code, atlas_str, atlas_ref, atlas_parcs, deriv_fold
                )
            elif proc_dict["method"] is None:
                self._process_method_none(supra, atlas_code, atlas_str, deriv_fold)
            elif proc_dict["method"] == "atlasbased":
                self._process_atlasbased(
                    supra, atlas_code, atlas_str, atlas_ref, deriv_fold, proc_dict, spam_thresh
                )

        self._accumulate_supra(supra)

    def _find_previous_parc(self, supra):
        """Locate a previously generated chimera parcellation for ``supra``."""
        chim_code = self.chim_code
        prev_parc = None
        supra_pos = list(self.parc_dict.keys()).index(supra)

        for i, prev_chim in enumerate(self.prev_chims):
            base_name = os.path.basename(prev_chim)
            ent_chim_dict = cltbids.str2entity(base_name)
            prev_chim_code = ent_chim_dict["atlas"].replace("chimera", "")

            if prev_chim_code[supra_pos] == chim_code[supra_pos]:
                prev_parc = prev_chim
                break

        return prev_parc

    def _load_supra_from_previous(self, supra, atlas_code, prev_parc):
        """Extract the supra-region parcellations from a previous chimera file."""
        # Get the position of the supra-region in the list of supra-regions
        prev_parc_obj = cltparc.Parcellation(parc_file=prev_parc)

        # Check if the supra-region has a left hemisphere parcellation
        if "lh" in self.supra_dict[supra][supra][atlas_code].keys():
            # Detecting the name to filter the parcellation for the left hemisphere
            lh_cad = self.supra_dict[supra][supra][atlas_code]["lh"]["name"][0]

            # Get the string until the second dash in the name. It will be used to
            # filter the parcellation for the left hemisphere
            lh_cad = "-".join(lh_cad.split("-")[:2])

            self.lh_supra_parc = copy.deepcopy(prev_parc_obj)
            self.lh_supra_parc.keep_by_name(names2keep=[lh_cad])

        # Check if the supra-region has a right hemisphere parcellation
        if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
            # Detecting the name to filter the parcellation for the right hemisphere
            rh_cad = self.supra_dict[supra][supra][atlas_code]["rh"]["name"][0]

            # Get the string until the second dash in the name.
            rh_cad = "-".join(rh_cad.split("-")[:2])

            self.rh_supra_parc = copy.deepcopy(prev_parc_obj)
            self.rh_supra_parc.keep_by_name(names2keep=[rh_cad])

        # Check if the supra-region has a mid hemisphere parcellation
        if "mid" in self.supra_dict[supra][supra][atlas_code].keys():
            # Detecting the name to filter the parcellation for the right hemisphere
            mid_cad = self.supra_dict[supra][supra][atlas_code]["mid"]["name"][0]

            # Get the string until the second dash in the name.
            mid_cad = "-".join(mid_cad.split("-")[:2])

            self.mid_supra_parc = copy.deepcopy(prev_parc_obj)
            self.mid_supra_parc.keep_by_name(names2keep=[mid_cad])

    def _process_comform2native(
        self, supra, atlas_code, atlas_str, atlas_ref, atlas_parcs, deriv_fold
    ):
        """Handle supra-regions whose method is ``comform2native``."""
        sub2proc = self.sub2proc
        deriv_dir = self.deriv_dir
        path_cad = self.path_cad
        fullid = self.fullid
        force = self.force
        cont_tech_freesurfer = self.cont_tech_freesurfer
        cont_image_freesurfer = self.cont_image_freesurfer
        atlas_src = self.parc_dict[supra]["source"]
        proc_dict = self.parc_dict[supra]["processing"]

        if proc_dict["labels"] == "freesurferextra":
            sub2proc.launch_freesurfer(
                force=force,
                extra_proc=supra.lower(),
                cont_tech=cont_tech_freesurfer,
                cont_image=cont_image_freesurfer,
            )

            fsextra_files = glob(
                os.path.join(
                    sub2proc.subjs_dir,
                    sub2proc.subj_id,
                    "mri",
                    "*" + atlas_parcs + ".mgz",
                )
            )
            if len(fsextra_files) == 0:
                raise ValueError("The Freesurfer extra parcellation was not found.")

            elif len(fsextra_files) > 1:
                lh_mgz_image = os.path.join(
                    sub2proc.subjs_dir,
                    sub2proc.subj_id,
                    "mri",
                    "lh." + atlas_parcs + ".mgz",
                )
                rh_mgz_image = os.path.join(
                    sub2proc.subjs_dir,
                    sub2proc.subj_id,
                    "mri",
                    "rh." + atlas_parcs + ".mgz",
                )
                lh_nii_image = os.path.join(
                    deriv_dir,
                    deriv_fold,
                    path_cad,
                    "anat",
                    fullid + "_hemi-L_atlas-" + atlas_str + "_dseg.nii.gz",
                )
                rh_nii_image = os.path.join(
                    deriv_dir,
                    deriv_fold,
                    path_cad,
                    "anat",
                    fullid + "_hemi-R_atlas-" + atlas_str + "_dseg.nii.gz",
                )

                if not os.path.isfile(lh_nii_image):
                    dir_name = os.path.dirname(lh_nii_image)
                    dir_name = Path(dir_name)
                    dir_name.mkdir(parents=True, exist_ok=True)

                    if atlas_ref == "conform":
                        sub2proc.conform2native(
                            mgz_conform=lh_mgz_image,
                            nii_native=lh_nii_image,
                            cont_image=cont_image_freesurfer,
                            cont_tech=cont_tech_freesurfer,
                            force=force,
                        )
                    else:
                        lh_nii_image = lh_mgz_image

                if not os.path.isfile(rh_nii_image):
                    dir_name = os.path.dirname(rh_nii_image)
                    dir_name = Path(dir_name)
                    dir_name.mkdir(parents=True, exist_ok=True)

                    if atlas_ref == "conform":
                        sub2proc.conform2native(
                            mgz_conform=rh_mgz_image,
                            nii_native=rh_nii_image,
                            cont_image=cont_image_freesurfer,
                            cont_tech=cont_tech_freesurfer,
                            force=force,
                        )
                    else:
                        rh_nii_image = rh_mgz_image

                lh_supra_parc = cltparc.Parcellation(parc_file=lh_nii_image)
                lh_supra_parc.index = self.supra_dict[supra][supra][atlas_code]["lh"][
                    "index"
                ]
                lh_supra_parc.name = self.supra_dict[supra][supra][atlas_code]["lh"][
                    "name"
                ]
                lh_supra_parc.color = self.supra_dict[supra][supra][atlas_code]["lh"][
                    "color"
                ]
                lh_supra_parc.keep_by_code(
                    codes2keep=self.supra_dict[supra][supra][atlas_code]["lh"]["index"]
                )
                lh_supra_parc.export_colortable(
                    out_file=lh_nii_image.replace(".nii.gz", ".lut"),
                    lut_type="lut",
                )
                lh_supra_parc.export_colortable(
                    out_file=lh_nii_image.replace(".nii.gz", ".tsv"),
                    lut_type="tsv",
                )
                self.lh_supra_parc = lh_supra_parc

                rh_supra_parc = cltparc.Parcellation(parc_file=rh_nii_image)
                rh_supra_parc.index = self.supra_dict[supra][supra][atlas_code]["rh"][
                    "index"
                ]
                rh_supra_parc.name = self.supra_dict[supra][supra][atlas_code]["rh"][
                    "name"
                ]
                rh_supra_parc.color = self.supra_dict[supra][supra][atlas_code]["rh"][
                    "color"
                ]
                rh_supra_parc.keep_by_code(
                    codes2keep=self.supra_dict[supra][supra][atlas_code]["rh"]["index"]
                )
                rh_supra_parc.export_colortable(
                    out_file=rh_nii_image.replace(".nii.gz", ".lut"),
                    lut_type="lut",
                )
                rh_supra_parc.export_colortable(
                    out_file=rh_nii_image.replace(".nii.gz", ".tsv"),
                    lut_type="tsv",
                )
                self.rh_supra_parc = rh_supra_parc

            elif len(fsextra_files) == 1:
                mgz_image = fsextra_files[0]
                nii_image = os.path.join(
                    deriv_dir,
                    deriv_fold,
                    path_cad,
                    "anat",
                    fullid + "_atlas-" + atlas_str + "_dseg.nii.gz",
                )
                if not os.path.isfile(nii_image):
                    dir_name = os.path.dirname(nii_image)
                    dir_name = Path(dir_name)
                    dir_name.mkdir(parents=True, exist_ok=True)

                    if atlas_ref == "conform":
                        sub2proc.conform2native(
                            mgz_conform=mgz_image,
                            nii_native=nii_image,
                            cont_image=cont_image_freesurfer,
                            cont_tech=cont_tech_freesurfer,
                            force=force,
                        )
                    else:
                        nii_image = mgz_image

                tmp_parc = cltparc.Parcellation(parc_file=nii_image)
                index, name, color, opacity = _mix_side_prop(
                    self.supra_dict[supra][supra][atlas_code]
                )
                tmp_parc.index = index
                tmp_parc.name = name
                tmp_parc.color = color
                tmp_parc.opacity = opacity

                tmp_parc.export_colortable(
                    out_file=nii_image.replace(".nii.gz", ".lut"),
                    lut_type="lut",
                )
                tmp_parc.export_colortable(
                    out_file=nii_image.replace(".nii.gz", ".tsv"),
                    lut_type="tsv",
                )

                # Left Hemisphere
                if "lh" in self.supra_dict[supra][supra][atlas_code].keys():
                    self.lh_supra_parc = copy.deepcopy(tmp_parc)
                    self.lh_supra_parc.keep_by_code(
                        codes2keep=self.supra_dict[supra][supra][atlas_code]["lh"][
                            "index"
                        ]
                    )

                # Right Hemisphere
                if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
                    self.rh_supra_parc = copy.deepcopy(tmp_parc)
                    self.rh_supra_parc.keep_by_code(
                        codes2keep=self.supra_dict[supra][supra][atlas_code]["rh"][
                            "index"
                        ]
                    )

                # Non-hemispheric structures
                if "mid" in self.supra_dict[supra][supra][atlas_code].keys():
                    self.mid_supra_parc = copy.deepcopy(tmp_parc)
                    self.mid_supra_parc.keep_by_code(
                        codes2keep=self.supra_dict[supra][supra][atlas_code]["mid"][
                            "index"
                        ]
                    )

        else:

            if atlas_src == "freesurfer":
                nii_image = os.path.join(
                    sub2proc.subjs_dir,
                    sub2proc.subj_id,
                    "tmp",
                    atlas_parcs + ".nii.gz",
                )
                mgz_image = os.path.join(
                    sub2proc.subjs_dir,
                    sub2proc.subj_id,
                    "mri",
                    atlas_parcs + ".mgz",
                )

            if self.aseg_parc is None:
                if not os.path.isfile(nii_image):
                    sub2proc.conform2native(
                        mgz_conform=mgz_image,
                        nii_native=nii_image,
                        cont_image=cont_image_freesurfer,
                        cont_tech=cont_tech_freesurfer,
                        force=force,
                    )
                    self.files2del.append(nii_image)

                tmp_parc = cltparc.Parcellation(parc_file=nii_image)

            else:
                tmp_parc = copy.deepcopy(self.aseg_parc)

            if supra == "Hypothalamus":
                hypo_parc = copy.deepcopy(tmp_parc)

                # Get the needed labels for the hypothalamus parcellation
                lh_vdc_index = self.supra_dict["AuxiliaryLUTTable"][supra][atlas_code][
                    "lh"
                ]["index"][0]
                rh_vdc_index = self.supra_dict["AuxiliaryLUTTable"][supra][atlas_code][
                    "rh"
                ]["index"][0]
                mid_third_vent = self.supra_dict["AuxiliaryLUTTable"]["Ventricle"][
                    atlas_code
                ]["mid"]["index"][0]

                # Keep only the third ventricle
                hypo_parc.keep_by_code(codes2keep=[mid_third_vent])

                morph = cltimg.MorphologicalOperations()
                dilated = morph.dilate_mm(
                    hypo_parc.data,
                    hypo_parc.affine,
                    shape="ball",
                    dilation_mm=5,
                )

                # Creating the vdc parcellation
                hypo_parc = copy.deepcopy(tmp_parc)
                hypo_parc.keep_by_code(codes2keep=[lh_vdc_index, rh_vdc_index])
                hypo_parc.data = dilated * hypo_parc.data
                tmp_parc_data = np.zeros_like(tmp_parc.data)

                tmp_parc_data[hypo_parc.data == lh_vdc_index] = 1
                tmp_parc_data[hypo_parc.data == rh_vdc_index] = 2

                tmp_parc = cltparc.Parcellation(
                    parc_file=tmp_parc_data,
                    affine=tmp_parc.affine,
                )

            # Left Hemisphere
            if "lh" in self.supra_dict[supra][supra][atlas_code].keys():
                lh_supra_parc = copy.deepcopy(tmp_parc)
                lh_supra_parc.index = self.supra_dict[supra][supra][atlas_code]["lh"][
                    "index"
                ]
                lh_supra_parc.name = self.supra_dict[supra][supra][atlas_code]["lh"][
                    "name"
                ]
                lh_supra_parc.color = self.supra_dict[supra][supra][atlas_code]["lh"][
                    "color"
                ]
                lh_supra_parc.keep_by_code(
                    codes2keep=self.supra_dict[supra][supra][atlas_code]["lh"]["index"]
                )
                self.lh_supra_parc = lh_supra_parc

            # Right Hemisphere
            if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
                rh_supra_parc = copy.deepcopy(tmp_parc)
                rh_supra_parc.index = self.supra_dict[supra][supra][atlas_code]["rh"][
                    "index"
                ]
                rh_supra_parc.name = self.supra_dict[supra][supra][atlas_code]["rh"][
                    "name"
                ]
                rh_supra_parc.color = self.supra_dict[supra][supra][atlas_code]["rh"][
                    "color"
                ]
                rh_supra_parc.keep_by_code(
                    codes2keep=self.supra_dict[supra][supra][atlas_code]["rh"]["index"]
                )
                self.rh_supra_parc = rh_supra_parc

            # Non-hemispheric structures
            if "mid" in self.supra_dict[supra][supra][atlas_code].keys():
                mid_supra_parc = copy.deepcopy(tmp_parc)
                mid_supra_parc.index = self.supra_dict[supra][supra][atlas_code]["mid"][
                    "index"
                ]
                mid_supra_parc.name = self.supra_dict[supra][supra][atlas_code]["mid"][
                    "name"
                ]
                mid_supra_parc.color = self.supra_dict[supra][supra][atlas_code]["mid"][
                    "color"
                ]
                mid_supra_parc.keep_by_code(
                    codes2keep=self.supra_dict[supra][supra][atlas_code]["mid"]["index"]
                )
                self.mid_supra_parc = mid_supra_parc

    def _process_method_none(self, supra, atlas_code, atlas_str, deriv_fold):
        """Handle supra-regions whose processing method is ``None`` (FSL FIRST)."""
        deriv_dir = self.deriv_dir
        path_cad = self.path_cad
        fullid = self.fullid
        t1 = self.t1
        force = self.force
        cont_tech_fsl = self.cont_tech_fsl
        cont_image_fsl = self.cont_image_fsl

        # Running FIRST if it is needed
        if atlas_code == "R":
            fsl_outdir = os.path.join(deriv_dir, deriv_fold, path_cad, "anat")
            first_nii = os.path.join(
                str(fsl_outdir), fullid + "_atlas-" + atlas_str + "_dseg.nii.gz"
            )

            if self.first_parc is None:
                fsl_outdir = Path(fsl_outdir)
                fsl_outdir.mkdir(parents=True, exist_ok=True)

                # Running the FIRST
                launch_fsl_first(
                    t1,
                    first_parc=first_nii,
                    cont_tech=cont_tech_fsl,
                    cont_image=cont_image_fsl,
                    force=force,
                )
                self.first_parc = cltparc.Parcellation(parc_file=first_nii)

            tmp_parc = copy.deepcopy(self.first_parc)
            index, name, color = _mix_side_prop(
                self.supra_dict[supra][supra][atlas_code]
            )
            tmp_parc.index = index
            tmp_parc.name = name
            tmp_parc.color = color

            tmp_parc.export_colortable(
                out_file=first_nii.replace(".nii.gz", ".lut"),
                lut_type="lut",
            )
            tmp_parc.export_colortable(
                out_file=first_nii.replace(".nii.gz", ".tsv"),
                lut_type="tsv",
            )
        # Left Hemisphere
        if "lh" in self.supra_dict[supra][supra][atlas_code].keys():
            self.lh_supra_parc = copy.deepcopy(tmp_parc)
            self.lh_supra_parc.keep_by_code(
                codes2keep=self.supra_dict[supra][supra][atlas_code]["lh"]["index"]
            )

        # Right Hemisphere
        if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
            self.rh_supra_parc = copy.deepcopy(tmp_parc)
            self.rh_supra_parc.keep_by_code(
                codes2keep=self.supra_dict[supra][supra][atlas_code]["rh"]["index"]
            )

        # Non-hemispheric structures
        if "mid" in self.supra_dict[supra][supra][atlas_code].keys():
            self.mid_supra_parc = copy.deepcopy(tmp_parc)
            self.mid_supra_parc.keep_by_code(
                codes2keep=self.supra_dict[supra][supra][atlas_code]["mid"]["index"]
            )

    def _process_atlasbased(
        self, supra, atlas_code, atlas_str, atlas_ref, deriv_fold, proc_dict, spam_thresh
    ):
        """Handle supra-regions built by registration to an atlas template."""
        deriv_dir = self.deriv_dir
        path_cad = self.path_cad
        fullid = self.fullid
        t1 = self.t1
        pipe_dict = self.pipe_dict
        cont_tech_ants = self.cont_tech_ants
        cont_image_ants = self.cont_image_ants
        ent_dict_fullid = self.ent_dict_fullid

        t1_temp = proc_dict["reference"]
        atlas = proc_dict["labels"]
        atlas_type = proc_dict["type"]

        # Basename for transformations
        spat_tf_ent = ent_dict_fullid.copy()
        spat_tf_ent = cltbids.delete_entity(
            spat_tf_ent, ["space", "desc", "suffix", "extension"]
        )
        spat_tf_ent["from"] = "T1w"
        spat_tf_ent["to"] = atlas_ref
        spat_tf_ent["mode"] = "image"
        spat_tf_ent["suffix"] = "xfm"
        spat_tf_ent["extension"] = "mat"
        xfm_base = os.path.join(
            deriv_dir,
            pipe_dict["outputs"]["transforms"],
            path_cad,
            "anat",
            cltbids.entity2str(spat_tf_ent),
        )
        work_dir = os.path.join(deriv_dir, deriv_fold, path_cad, "anat")

        # Create the working directory if it does not exist using the pathlib library
        work_dir = Path(work_dir)

        out_parc_spam = os.path.join(
            str(work_dir), fullid + "_atlas-" + atlas_str + "_probseg.nii.gz"
        )
        out_parc_maxp = os.path.join(
            str(work_dir), fullid + "_atlas-" + atlas_str + "_dseg.nii.gz"
        )
        out_parc_lut = out_parc_maxp.replace(".nii.gz", ".lut")
        out_parc_tsv = out_parc_maxp.replace(".nii.gz", ".tsv")

        if (
            not os.path.isfile(out_parc_maxp)
            or not os.path.isfile(out_parc_lut)
            or not os.path.isfile(out_parc_tsv)
            or self.force
        ):
            work_dir.mkdir(parents=True, exist_ok=True)

            # Detecting the side
            sides_ids = list(self.supra_dict[supra][supra][atlas_code].keys())
            sides_ids = sorted(
                sides_ids, key=lambda x: not ("lh" in x or "rh" in x)
            )

            # Masking the cerebellum from T1w image
            tmp_t1 = t1
            if supra == "Cerebellum":
                if self.parc_dict[supra]["name"] == "SUIT":
                    tmp_t1 = os.path.join(str(work_dir), "tmp_cerb_bs.nii.gz")
                    self.files2del.append(tmp_t1)

                    # Hard coding the cerebellum labels for the aseg parcellation
                    cltimg.crop_image_from_mask(
                        t1,
                        self.aseg_parc.data,
                        tmp_t1,
                        [7, 8, 16, 46, 47, 15, 16],
                    )

            if atlas_type == "spam":
                cltseg.abased_parcellation(
                    tmp_t1,
                    t1_temp,
                    atlas,
                    out_parc_spam,
                    xfm_base,
                    cont_tech=cont_tech_ants,
                    cont_image=cont_image_ants,
                )

                if supra == "Cerebellum":
                    if self.parc_dict[supra]["name"] == "SUIT":
                        os.remove(tmp_t1)
                        cltimg.cropped_to_native(out_parc_spam, t1, out_parc_spam)

                for side_cont, side in enumerate(sides_ids):
                    vol_indexes = (
                        np.array(
                            self.supra_dict[supra][supra][atlas_code][side]["index"]
                        )
                        - 1
                    )
                    tmp_par_file = os.path.join(
                        str(work_dir),
                        fullid
                        + "_hemi-"
                        + side
                        + "_atlas-"
                        + atlas_str
                        + "_dseg.nii.gz",
                    )
                    self.files2del.append(tmp_par_file)

                    tmp_parc_file = cltimg.spams2maxprob(
                        out_parc_spam,
                        prob_thresh=spam_thresh,
                        vol_indexes=vol_indexes,
                        maxp_name=tmp_par_file,
                    )

                    tmp_parc = cltparc.Parcellation(parc_file=tmp_parc_file)
                    tmp_parc.index = vol_indexes + 1
                    tmp_parc.name = self.supra_dict[supra][supra][atlas_code][side][
                        "name"
                    ]
                    tmp_parc.color = self.supra_dict[supra][supra][atlas_code][side][
                        "color"
                    ]

                    if side in self.supra_dict[supra][supra]["F"].keys():
                        aseg_code = self.supra_dict[supra][supra]["F"][side]["index"]
                        tmp_parc.apply_mask(
                            image_mask=self.aseg_parc,
                            mask_codes=aseg_code,
                            fill=True,
                        )

                    else:
                        # Selecting all the region in case there is no definition of
                        # left and right hemispheres
                        all_reg_codes = []
                        for side_f in self.supra_dict[supra][supra]["F"].keys():
                            aseg_code = self.supra_dict[supra][supra]["F"][side_f][
                                "index"
                            ]
                            all_reg_codes = all_reg_codes + aseg_code

                        glob_mask_parc = copy.deepcopy(self.aseg_parc)
                        group_dict = {
                            1: {
                                "index": all_reg_codes,
                                "name": "all_aseg",
                                "color": "#FFFFFF",
                                "opacity": 1,
                            }
                        }
                        glob_mask_parc.group_by_codes(group_dict)
                        tmp_parc.apply_mask(image_mask=glob_mask_parc, mask_codes=1)

                    # Adjusting the values to the ones existing on the 3D image
                    tmp_parc.adjust_values()

                    if side_cont == 0:
                        def_parc = copy.deepcopy(tmp_parc)
                    else:
                        def_parc.add_parcellation(tmp_parc)

                    # Removing the temporal side images
                    if os.path.isfile(tmp_parc_file):
                        os.remove(tmp_parc_file)

                def_parc.save_parcellation(
                    out_file=out_parc_maxp,
                    affine=def_parc.affine,
                    lut_type=["lut", "tsv"],
                )

            elif atlas_type == "maxprob":

                cltseg.abased_parcellation(
                    tmp_t1,
                    t1_temp,
                    atlas,
                    out_parc_maxp,
                    xfm_base,
                    atlas_type="maxprob",
                    cont_tech=cont_tech_ants,
                    cont_image=cont_image_ants,
                )

                if supra == "Cerebellum":
                    if self.parc_dict[supra]["name"] == "SUIT":
                        os.remove(tmp_t1)
                        cltimg.cropped_to_native(out_parc_maxp, t1, out_parc_maxp)

                for side_cont, side in enumerate(sides_ids):
                    tmp_par_file = os.path.join(
                        work_dir,
                        fullid
                        + "_hemi-"
                        + side
                        + "_atlas-"
                        + atlas_str
                        + "_dseg.nii.gz",
                    )
                    self.files2del.append(tmp_par_file)

                    tmp_parc = cltparc.Parcellation(parc_file=out_parc_maxp)
                    tmp_parc.data = np.round(tmp_parc.data)
                    tmp_parc.index = np.array(
                        self.supra_dict[supra][supra][atlas_code][side]["index"]
                    )
                    tmp_parc.name = self.supra_dict[supra][supra][atlas_code][side][
                        "name"
                    ]
                    tmp_parc.color = self.supra_dict[supra][supra][atlas_code][side][
                        "color"
                    ]
                    tmp_parc.opacity = self.supra_dict[supra][supra][atlas_code][side][
                        "opacity"
                    ]

                    if side in self.supra_dict[supra][supra]["F"].keys():
                        aseg_code = self.supra_dict[supra][supra]["F"][side]["index"]
                        tmp_parc.apply_mask(
                            image_mask=self.aseg_parc,
                            mask_codes=aseg_code,
                            fill=True,
                        )

                    else:
                        # Selecting all the region in case there is no definition of
                        # left and right hemispheres
                        all_reg_codes = []
                        for side_f in self.supra_dict[supra][supra]["F"].keys():
                            aseg_code = self.supra_dict[supra][supra]["F"][side_f][
                                "index"
                            ]
                            all_reg_codes = all_reg_codes + aseg_code

                        glob_mask_parc = copy.deepcopy(self.aseg_parc)
                        glob_mask_parc.group_by_code(
                            codes2group=all_reg_codes, new_codes=1
                        )
                        tmp_parc.apply_mask(image_mask=glob_mask_parc, mask_codes=1)

                    # Adjusting the values to the ones existing on the 3D image
                    tmp_parc.adjust_values()
                    if side_cont == 0:
                        def_parc = copy.deepcopy(tmp_parc)
                    else:
                        def_parc.add_parcellation(tmp_parc)

                def_parc.save_parcellation(
                    out_file=out_parc_maxp,
                    affine=def_parc.affine,
                    lut_type=["lut", "tsv"],
                )

        tmp_parc = cltparc.Parcellation(parc_file=out_parc_maxp)
        index, name, color, opacity = _mix_side_prop(
            self.supra_dict[supra][supra][atlas_code]
        )
        tmp_parc.index = index
        tmp_parc.name = name
        tmp_parc.color = color
        tmp_parc.opacity = opacity
        tmp_parc.adjust_values()

        tmp_parc.export_colortable(
            out_file=out_parc_maxp.replace(".nii.gz", ".lut"),
            lut_type="lut",
        )
        tmp_parc.export_colortable(
            out_file=out_parc_maxp.replace(".nii.gz", ".tsv"),
            lut_type="tsv",
        )
        # Left Hemisphere
        if "lh" in self.supra_dict[supra][supra][atlas_code].keys():
            self.lh_supra_parc = copy.deepcopy(tmp_parc)
            self.lh_supra_parc.keep_by_code(
                codes2keep=self.supra_dict[supra][supra][atlas_code]["lh"]["index"]
            )

        # Right Hemisphere
        if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
            self.rh_supra_parc = copy.deepcopy(tmp_parc)
            self.rh_supra_parc.keep_by_code(
                codes2keep=self.supra_dict[supra][supra][atlas_code]["rh"]["index"]
            )

        # Non-hemispheric structures
        if "mid" in self.supra_dict[supra][supra][atlas_code].keys():
            self.mid_supra_parc = copy.deepcopy(tmp_parc)
            self.mid_supra_parc.keep_by_code(
                codes2keep=self.supra_dict[supra][supra][atlas_code]["mid"]["index"]
            )

    def _accumulate_supra(self, supra):
        """Fold the supra-region parcellations into the global buffers."""
        if self.lh_supra_parc is not None:
            self.lh_supra_parc.rearrange()

            if "F" in self.supra_dict[supra][supra].keys():
                # Use the FreeSurfer parcellation to detect the voxels that are not
                # in the lh_supra_parc
                lh2refill_parc = copy.deepcopy(self.aseg_parc)
                lh2refill_parc.keep_by_code(
                    codes2keep=self.supra_dict[supra][supra]["F"]["lh"]["index"]
                )

                # Find the voxels that are not in the lh_supra_parc and are in the lh2refill
                ind = np.where(
                    (self.lh_supra_parc.data == 0) & (lh2refill_parc.data != 0)
                )
                self.lh2refill[ind] = 1
                del lh2refill_parc

            # Add the parcellation to the global left subcortical parcellation
            self.lh_parc.add_parcellation(self.lh_supra_parc, append=True)
            self.nlh_subc = len(self.lh_parc.index)
            self.lh_supra_parc = None

        if self.rh_supra_parc is not None:
            self.rh_supra_parc.rearrange()

            if "F" in self.supra_dict[supra][supra].keys():
                # Use the FreeSurfer parcellation to detect the voxels that are not
                # in the rh_supra_parc
                rh2refill_parc = copy.deepcopy(self.aseg_parc)
                rh2refill_parc.keep_by_code(
                    codes2keep=self.supra_dict[supra][supra]["F"]["rh"]["index"]
                )

                # Find the voxels that are not in the rh_supra_parc and are in the rh2refill
                ind = np.where(
                    (self.rh_supra_parc.data == 0) & (rh2refill_parc.data != 0)
                )
                self.rh2refill[ind] = 1
                del rh2refill_parc

            # Add the parcellation to the global right subcortical parcellation
            self.rh_parc.add_parcellation(self.rh_supra_parc, append=True)
            self.nrh_subc = len(self.rh_parc.index)
            self.rh_supra_parc = None

        if self.mid_supra_parc is not None:
            self.mid_supra_parc.rearrange()
            self.mid_parc.add_parcellation(self.mid_supra_parc, append=True)
            self.mid_supra_parc = None

            if self.rh2refill is not None:
                #  Remove the voxels that are in the mid_parc and are in the rh2refill
                ind = np.where((self.mid_parc.data != 0) & (self.rh2refill != 0))
                self.rh2refill[ind] = 0

            if self.lh2refill is not None:
                #  Remove the voxels that are in the mid_parc and are in the lh2refill
                ind = np.where((self.mid_parc.data != 0) & (self.lh2refill != 0))
                self.lh2refill[ind] = 0

    def _count_subcortical_regions(self):
        """Count the regions accumulated in each hemisphere buffer."""
        self.nlh_subc = len(self.lh_parc.index)
        self.nrh_subc = len(self.rh_parc.index)
        self.nmid_subc = len(self.mid_parc.index)

    # ------------------------------------------------- step: final assembly
    def _assemble_and_write(self):
        """Assemble the final parcellation image(s) and write them to disk."""
        date_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        if self.bool_ctx:
            self._assemble_with_cortex(date_time)
        else:
            self._assemble_without_cortex(date_time)

    def _assemble_with_cortex(self, date_time):
        """Assemble and write parcellations that include a cortical atlas."""
        deriv_dir = self.deriv_dir
        path_cad = self.path_cad
        fullid = self.fullid
        chim_code = self.chim_code
        chim_dir = self.out_dir
        bool_mixwm = self.bool_mixwm
        force = self.force
        sub2proc = self.sub2proc
        prev_chims = self.prev_chims
        cont_tech_freesurfer = self.cont_tech_freesurfer
        cont_image_freesurfer = self.cont_image_freesurfer
        t1_image = self.t1_image
        affine = self.affine

        # Atributes for the cortical parcellation
        atlas_names = self.parc_dict["Cortical"]["parcels"]

        proc_dict = self.parc_dict["Cortical"]["processing"]
        ctx_meth = proc_dict["method"]

        nctx_parc = len(self.parc_dict["Cortical"]["processing"]["labels"]["lh"])
        for c in np.arange(nctx_parc):

            # Temporal header lines
            glob_header_info_tmp = copy.deepcopy(self.glob_header_info)

            ## -------- Cortical parcellation for the left hemisphere -----------
            # Creating the name for the output file
            lh_in_parc = self.parc_dict["Cortical"]["processing"]["labels"]["lh"][c]
            at_name = [s for s in atlas_names if s in lh_in_parc]
            lh_out_annot = os.path.join(
                deriv_dir,
                self.parc_dict["Cortical"]["deriv_surffold"],
                path_cad,
                "anat",
                fullid + "_hemi-L" + "_" + "".join(at_name) + "_dseg.annot",
            )

            ## -------- Cortical parcellation for the right hemisphere ----------
            # Creating the name for the output file
            rh_in_parc = self.parc_dict["Cortical"]["processing"]["labels"]["rh"][c]
            at_name = [s for s in atlas_names if s in rh_in_parc]
            at_name = "".join(at_name)
            rh_out_annot = os.path.join(
                deriv_dir,
                self.parc_dict["Cortical"]["deriv_surffold"],
                path_cad,
                "anat",
                fullid + "_hemi-R" + "_" + at_name + "_dseg.annot",
            )

            if ctx_meth == "annot2indiv":
                # Moving to individual space
                sub2proc.annot2ind(
                    ref_id=self.parc_dict["Cortical"]["processing"]["reference"],
                    hemi="lh",
                    fs_annot=lh_in_parc,
                    ind_annot=lh_out_annot,
                    cont_tech=cont_tech_freesurfer,
                    cont_image=cont_image_freesurfer,
                    force=force,
                )

                sub2proc.annot2ind(
                    ref_id=self.parc_dict["Cortical"]["processing"]["reference"],
                    hemi="rh",
                    fs_annot=rh_in_parc,
                    ind_annot=rh_out_annot,
                    cont_tech=cont_tech_freesurfer,
                    cont_image=cont_image_freesurfer,
                    force=force,
                )

            if ctx_meth == "gcs2indiv":
                # Moving to individual space
                sub2proc.gcs2ind(
                    fs_gcs=lh_in_parc,
                    hemi="lh",
                    ind_annot=lh_out_annot,
                    cont_tech=cont_tech_freesurfer,
                    cont_image=cont_image_freesurfer,
                    force=force,
                )

                sub2proc.gcs2ind(
                    fs_gcs=rh_in_parc,
                    hemi="rh",
                    ind_annot=rh_out_annot,
                    cont_tech=cont_tech_freesurfer,
                    cont_image=cont_image_freesurfer,
                    force=force,
                )

            # Copying to the labels folder
            temp_lh = os.path.join(
                sub2proc.subjs_dir,
                sub2proc.subj_id,
                "label",
                "lh." + at_name + ".annot",
            )
            shutil.copyfile(lh_out_annot, temp_lh)

            # Copying to the labels folder
            temp_rh = os.path.join(
                sub2proc.subjs_dir,
                sub2proc.subj_id,
                "label",
                "rh." + at_name + ".annot",
            )
            shutil.copyfile(rh_out_annot, temp_rh)

            ## -------- Creating the volumetric parcellation --------------------
            out_vol_dir = os.path.join(
                deriv_dir,
                self.parc_dict["Cortical"]["deriv_volfold"],
                path_cad,
                "anat",
            )
            if self.growwm is None:
                self.growwm = ["0"]

            ent_dict = cltbids.str2entity(at_name)
            if "scale" in ent_dict.keys():
                scale_cad = "Scale: {}".format(ent_dict["scale"])
            else:
                scale_cad = None

            if "seg" in ent_dict.keys():
                seg_cad = "Segmentation: {}".format(ent_dict["seg"])
            else:
                seg_cad = None

            if scale_cad is not None or seg_cad is not None:
                if scale_cad is not None and seg_cad is not None:
                    cad2add = ". " + scale_cad + " - " + seg_cad

                elif scale_cad is not None and seg_cad is None:
                    cad2add = ". " + scale_cad

                elif scale_cad is None and seg_cad is not None:
                    cad2add = ". " + seg_cad

                glob_header_info_tmp[0] = glob_header_info_tmp[0] + cad2add

            for ngrow in np.arange(len(self.growwm)):
                if self.growwm[ngrow] == "0":
                    out_vol_name = fullid + "_" + at_name + "_dseg.nii.gz"
                else:

                    ent_dict = cltbids.str2entity(at_name)
                    if "desc" in ent_dict.keys():
                        if bool_mixwm:
                            ent_dict["desc"] = (
                                ent_dict["desc"]
                                + "grow"
                                + str(self.growwm[ngrow])
                                + "mm+mixwm"
                            )
                        else:
                            ent_dict["desc"] = (
                                ent_dict["desc"]
                                + "grow"
                                + str(self.growwm[ngrow])
                                + "mm"
                            )
                        tmp_str = cltbids.entity2str(ent_dict)
                        out_vol_name = fullid + "_" + tmp_str + "_dseg.nii.gz"
                    else:
                        if bool_mixwm:
                            out_vol_name = (
                                fullid
                                + "_"
                                + at_name
                                + "_desc-grow"
                                + str(self.growwm[ngrow])
                                + "mm+mixwm_dseg.nii.gz"
                            )
                        else:
                            out_vol_name = (
                                fullid
                                + "_"
                                + at_name
                                + "_desc-grow"
                                + str(self.growwm[ngrow])
                                + "mm_dseg.nii.gz"
                            )

                # Get the final name for the cortical parcellation
                chim_parc_name = cltbids.replace_entity_value(
                    out_vol_name, {"atlas": "chimera" + chim_code}
                )
                chim_parc_file = os.path.join(str(chim_dir), chim_parc_name)
                chim_parc_lut = os.path.join(
                    str(chim_dir),
                    cltbids.replace_entity_value(
                        chim_parc_name, {"extension": "lut"}
                    ),
                )
                chim_parc_tsv = os.path.join(
                    str(chim_dir),
                    cltbids.replace_entity_value(
                        chim_parc_name, {"extension": "tsv"}
                    ),
                )
                if (
                    not os.path.isfile(chim_parc_file)
                    or not os.path.isfile(chim_parc_lut)
                    or not os.path.isfile(chim_parc_tsv)
                    or force
                ):
                    self._compose_and_save_cortical(
                        c,
                        ngrow,
                        at_name,
                        out_vol_name,
                        out_vol_dir,
                        chim_parc_name,
                        chim_parc_file,
                        chim_parc_lut,
                        glob_header_info_tmp,
                        date_time,
                        prev_chims,
                    )

    def _compose_and_save_cortical(
        self,
        c,
        ngrow,
        at_name,
        out_vol_name,
        out_vol_dir,
        chim_parc_name,
        chim_parc_file,
        chim_parc_lut,
        glob_header_info_tmp,
        date_time,
        prev_chims,
    ):
        """Compose and save a single cortical Chimera parcellation image."""
        chim_code = self.chim_code
        bool_mixwm = self.bool_mixwm
        sub2proc = self.sub2proc
        cont_tech_freesurfer = self.cont_tech_freesurfer
        cont_image_freesurfer = self.cont_image_freesurfer
        t1_image = self.t1_image
        affine = self.affine

        # Creating the first part of the headers
        part_header = [f"# $Id: {chim_parc_lut} {date_time} \n"]

        part_header.append(
            f"# Corresponding parcellation: {chim_parc_file} \n"
        )

        lut_header = part_header + glob_header_info_tmp

        ent_tmp_dict = cltbids.str2entity(cltmisc.get_real_basename(chim_parc_name))
        # Check if there is a previous cortical parcellation

        ctx_prev_chims = cltmisc.filter_by_substring(
            prev_chims, "chimera" + chim_code[0]
        )
        if len(ctx_prev_chims) > 0:
            if "scale" in ent_tmp_dict.keys():
                ctx_prev_chims = cltmisc.filter_by_substring(
                    ctx_prev_chims, "scale-" + ent_tmp_dict["scale"]
                )
        if len(ctx_prev_chims) > 0:
            if "seg" in ent_tmp_dict.keys():
                ctx_prev_chims = cltmisc.filter_by_substring(
                    ctx_prev_chims, "seg-" + ent_tmp_dict["seg"]
                )

        prev_parc = None
        if len(ctx_prev_chims) > 0:
            for i, prev_chim in enumerate(ctx_prev_chims):
                base_name = cltmisc.get_real_basename(prev_chim)
                ent_chim_dict = cltbids.str2entity(base_name)

                if (
                    "desc" in ent_tmp_dict.keys()
                ):  # If desc is in the name of the desired parcellation we will look for a previous parcellation with the same desc
                    if (
                        "desc" in ent_chim_dict.keys()
                        and ent_chim_dict["desc"] == ent_tmp_dict["desc"]
                    ):
                        prev_parc = prev_chim
                        prev_tmp_parc = cltparc.Parcellation(parc_file=prev_parc)
                        break
                elif (
                    "desc" not in ent_tmp_dict.keys()
                ):  # If desc is not in the name of the desired parcellation
                    if (
                        "desc" not in ent_chim_dict.keys()
                    ):  # We will look for a previous parcellation without desc
                        prev_parc = prev_chim
                        prev_tmp_parc = cltparc.Parcellation(parc_file=prev_parc)
                        break
                    else:  # Or if it has desc
                        if (
                            "mixwm" not in ent_chim_dict["desc"]
                        ):  # but it is not a mixed white matter parcellation we will also consider it as a previous parcellation
                            prev_parc = prev_chim
                            prev_tmp_parc = cltparc.Parcellation(parc_file=prev_parc)

                            # but we will change the values of the cortical white
                            # matter to 3000 to match the values of the current
                            # parcellation
                            ind_ctx_wm = np.where(
                                (prev_tmp_parc.data > 3000)
                                & (prev_tmp_parc.data < 4000)
                            )
                            prev_tmp_parc.data[ind_ctx_wm] = 3000
                            prev_tmp_parc.adjust_values()
                            break

        # Creating the joined parcellation
        ref_image = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
        chim_parc = cltparc.Parcellation(parc_file=ref_image, affine=affine)

        if prev_parc is not None:
            ctx_parc = copy.deepcopy(prev_tmp_parc)

            brain_wm_parc = copy.deepcopy(prev_tmp_parc)
            brain_wm_parc.keep_by_name(names2keep=["wm-brain-whitematter"])

        else:

            sub2proc.surf2vol(
                atlas=at_name,
                out_vol=os.path.join(out_vol_dir, out_vol_name),
                gm_grow=self.growwm[ngrow],
                bool_mixwm=bool_mixwm,
                force=False,
                bool_native=True,
                color_table=["tsv", "lut"],
                cont_tech=cont_tech_freesurfer,
                cont_image=cont_image_freesurfer,
            )

            ctx_parc = cltparc.Parcellation(
                parc_file=os.path.join(out_vol_dir, out_vol_name)
            )
            ctx_parc.remove_by_name(
                names2remove=[
                    "unknown",
                    "medialwall",
                    "corpuscallosum",
                ]
            )

            # Detect the global White Matter
            brain_wm_parc = copy.deepcopy(ctx_parc)
            brain_wm_parc.keep_by_code(codes2keep=[2, 41, 5001, 5002, 7, 46])
            ind = np.where(brain_wm_parc.data != 0)
            brain_wm_parc.data[ind] = 1
            brain_wm_parc.index = [1]
            brain_wm_parc.name = ["wm-brain-whitematter"]
            brain_wm_parc.color = ["#ffffff"]
            brain_wm_parc.opacity = [1.0]
            brain_wm_parc.rearrange(offset=2999)
            brain_wm_parc.data[np.where(self.lh2refill)] = 3000
            brain_wm_parc.data[np.where(self.rh2refill)] = 3000

        lh_ctx_parc = copy.deepcopy(ctx_parc)
        rh_ctx_parc = copy.deepcopy(ctx_parc)

        lh_ctx_parc.keep_by_name(names2keep="ctx-lh-")
        nlh_ctx = len(lh_ctx_parc.index)
        rh_ctx_parc.keep_by_name(names2keep="ctx-rh-")
        nrh_ctx = len(rh_ctx_parc.index)

        # White Matter for the Right Hemisphere
        tmp_rh = cltmisc.filter_by_substring(ctx_parc.name, "wm-rh-")
        if tmp_rh:
            rh_wm_parc = copy.deepcopy(ctx_parc)
            rh_wm_parc.keep_by_name(names2keep=tmp_rh)
            rh_wm_parc.rearrange(offset=3000)

        # White Matter for the Left Hemisphere
        tmp_lh = cltmisc.filter_by_substring(ctx_parc.name, "wm-lh-")
        if tmp_lh:
            lh_wm_parc = copy.deepcopy(ctx_parc)
            lh_wm_parc.keep_by_name(names2keep=tmp_lh)
            lh_wm_parc.rearrange(offset=3000 + nrh_ctx + self.nrh_subc)

        # Adding the right cortical parcellation to the final image
        rh_ctx_parc.rearrange()
        chim_parc.add_parcellation(rh_ctx_parc, append=True)
        del rh_ctx_parc

        # Adding the right non-cortical parcellation to the final image
        if self.rh_parc is not None:
            chim_parc.add_parcellation(self.rh_parc, append=True)

        # Adding the left cortical parcellation to the final image
        lh_ctx_parc.rearrange()
        chim_parc.add_parcellation(lh_ctx_parc, append=True)
        del lh_ctx_parc

        # Adding the left non-cortical parcellation to the final image
        if self.lh_parc is not None:
            chim_parc.add_parcellation(self.lh_parc, append=True)

        # Adding the regions that do not belong to any hemisphere to the final image
        if self.mid_parc is not None:
            chim_parc.add_parcellation(self.mid_parc, append=True)

        # Adding the white matter to the final image
        chim_parc.add_parcellation(brain_wm_parc, append=False)
        del brain_wm_parc

        if tmp_rh:
            chim_parc.add_parcellation(rh_wm_parc, append=False)
            del rh_wm_parc

        if tmp_lh:
            chim_parc.add_parcellation(lh_wm_parc, append=False)
            del lh_wm_parc

        # Adding the extra regions
        if self.extra_parc is not None:
            # Detecting if there is region overlap and removing it
            tmp_extra = copy.deepcopy(self.extra_parc)
            mask = np.logical_and(tmp_extra.data != 0, chim_parc.data != 0)
            indexes = np.where(mask)
            tmp_extra.data[indexes] = 0
            tmp_extra.adjust_values()
            chim_parc.add_parcellation(tmp_extra, append=False)
            del tmp_extra

        # Saving the FINAL parcellation
        chim_parc.save_parcellation(
            out_file=chim_parc_file,
            affine=affine,
            headerlines=lut_header,
            lut_type=["lut", "tsv"],
        )
        del chim_parc

    def _assemble_without_cortex(self, date_time):
        """Assemble and write the parcellation when no cortical atlas is used."""
        fullid = self.fullid
        chim_code = self.chim_code
        chim_dir = self.out_dir
        force = self.force
        t1_image = self.t1_image
        affine = self.affine

        out_vol_name = fullid + "_dseg.nii.gz"

        chim_parc_name = cltbids.insert_entity(
            out_vol_name, {"atlas": "chimera" + chim_code}
        )
        chim_parc_file = os.path.join(str(chim_dir), chim_parc_name)
        chim_parc_lut = os.path.join(
            str(chim_dir),
            cltbids.replace_entity_value(chim_parc_name, {"extension": "lut"}),
        )
        chim_parc_tsv = os.path.join(
            str(chim_dir),
            cltbids.replace_entity_value(chim_parc_name, {"extension": "tsv"}),
        )

        if (
            not os.path.isfile(chim_parc_file)
            or not os.path.isfile(chim_parc_lut)
            or not os.path.isfile(chim_parc_tsv)
            or force
        ):
            part_header = [f"# $Id: {chim_parc_lut} {date_time} \n"]
            part_header.append(
                f"# Corresponding parcellation: {chim_parc_file} \n"
            )

            lut_header = part_header + self.glob_header_info
            lut_header = lut_header + ["\n"]

            # Creating the joined parcellation
            ref_image = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
            chim_parc = cltparc.Parcellation(parc_file=ref_image, affine=affine)

            # Adding the right non-cortical parcellation to the final image
            if self.rh_parc is not None:
                self.rh_parc.rearrange()
                chim_parc.add_parcellation(self.rh_parc, append=True)

            # Adding the left non-cortical parcellation to the final image
            if self.lh_parc is not None:
                self.lh_parc.rearrange()
                chim_parc.add_parcellation(self.lh_parc, append=True)

            # Adding the regions that do not belong to any hemisphere to the final image
            if self.mid_parc is not None:
                self.mid_parc.rearrange()
                chim_parc.add_parcellation(self.mid_parc, append=True)

            # Saving the FINAL parcellation
            chim_parc.save_parcellation(
                out_file=chim_parc_file,
                affine=affine,
                headerlines=lut_header,
                lut_type=["lut", "tsv"],
            )
            del chim_parc

    # ------------------------------------------------------- step: cleanup
    def _cleanup(self):
        """Remove the intermediate/temporary files collected during processing."""
        for tmp_file in self.files2del:
            if os.path.isfile(tmp_file):
                os.remove(tmp_file)
