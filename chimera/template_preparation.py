#!/usr/bin/env python3
"""Template preparation for the Chimera parcellation.

This module decomposes the original ``Chimera.prepare_templates`` monolith into
small, named functions.  The behaviour is preserved exactly; the only structural
change is that ``pipe_dict`` is threaded explicitly instead of being read from a
module-level global.

The public entry point is :func:`prepare_templates`, which mutates
``chim.parc_dict[supra]["processing"]`` for every supra-region, exactly like the
original method did.
"""

import os
from glob import glob
from os import PathLike

import clabtoolkit.freesurfertools as cltfree
import clabtoolkit.misctools as cltmisc
from templateflow import api as tflow

from .config_manager import _pipeline_info, _set_templateflow_home


def prepare_templates(chim, fssubj_dir: str = None, pipe_dict: dict = None):
    """Prepare the templates for the Chimera parcellation.

    Based on the parcellation code, this downloads the necessary templates from
    the TemplateFlow repository or sets up the templates directory for each
    supra-region.  It also creates the necessary symlinks to the FreeSurfer
    directory if the cortical parcellation is selected.

    Parameters
    ----------
    chim : Chimera
        The Chimera object whose ``parc_dict`` will be populated with a
        ``processing`` entry per supra-region.

    fssubj_dir : str
        FreeSurfer directory.

    pipe_dict : dict
        Pipeline configuration dictionary.  If ``None`` it is resolved with
        :func:`_pipeline_info`.
    """

    if pipe_dict is None:
        pipe_dict = _pipeline_info()

    # Setting up the FreeSurfer directory
    cltfree.FreeSurferSubject.set_freesurfer_directory(fssubj_dir)

    # Create the simlink to the FreeSurfer directory
    if "Cortical" in chim.parc_dict.keys():
        cltfree.create_fsaverage_links(
            fssubj_dir,
            fsavg_dir=None,
            refsubj_name=chim.parc_dict["Cortical"]["reference"],
        )

    # Detecting the base directory
    chim_dir = os.path.dirname(os.path.abspath(__file__))

    # Reading the names of the supra-regions
    supra_names = list(chim.parc_dict.keys())

    for supra in supra_names:
        if supra == "Cortical":
            meth_dict = _prepare_cortical_templates(
                chim, supra, fssubj_dir, chim_dir, pipe_dict
            )
        else:
            meth_dict = _prepare_volumetric_templates(chim, supra, chim_dir, pipe_dict)

        chim.parc_dict[supra]["processing"] = meth_dict


def _prepare_cortical_templates(chim, supra, fssubj_dir, chim_dir, pipe_dict):
    """Prepare the cortical surface templates for a supra-region."""

    atlas_src = chim.parc_dict[supra]["source"]
    atlas_str = chim.parc_dict[supra]["atlas"]
    atlas_ref = chim.parc_dict[supra]["reference"]

    # Atributes for the cortical parcellation
    atlas_type = chim.parc_dict[supra]["type"]
    atlas_names = chim.parc_dict[supra]["parcels"]

    # Selecting the source and downloading the parcellation
    if atlas_src == "templateflow":
        atlas_ext = ".gii"
        method = "annot2indiv"

        tflow_home = _set_templateflow_home(
            pipe_dict["packages"]["templateflow"]["home_dir"]
        )
        ctx_parc_lh = tflow.get(
            template=atlas_ref,
            atlas=atlas_str,
            hemi="L",
            suffix="dseg",
            extension=".label.gii",
        )
        ctx_parc_rh = tflow.get(
            template=atlas_ref,
            atlas=atlas_str,
            hemi="R",
            suffix="dseg",
            extension=".label.gii",
        )

        # Convert the list of PosfixPath to strings
        if isinstance(ctx_parc_lh, list):
            ctx_parc_lh = [str(x) for x in ctx_parc_lh]

        elif isinstance(ctx_parc_lh, PathLike):
            ctx_parc_lh = [str(ctx_parc_lh)]

        if isinstance(ctx_parc_rh, list):
            ctx_parc_rh = [str(x) for x in ctx_parc_rh]

        elif isinstance(ctx_parc_rh, PathLike):
            ctx_parc_rh = [str(ctx_parc_rh)]

        # Select the files that contain the atlas names
        ctx_parc_lh = cltmisc.filter_by_substring(
            ctx_parc_lh, atlas_names, bool_case=False
        )
        ctx_parc_rh = cltmisc.filter_by_substring(
            ctx_parc_rh, atlas_names, bool_case=False
        )
        ctx_parc_lh.sort()
        ctx_parc_rh.sort()

        ctx_parc_lh_annot = []
        ctx_parc_rh_annot = []
        for i, parc_file in enumerate(ctx_parc_lh):

            # Detect which element in atlas_names is in the string ctx_parc_lh
            at_name = [s for s in atlas_names if s in ctx_parc_lh[i]]

            if at_name:

                # Moving the gifti to native space
                tmp_annot = os.path.join(
                    fssubj_dir,
                    atlas_ref,
                    "label",
                    "lh." + at_name[0] + ".annot",
                )
                tmp_refsurf = os.path.join(
                    fssubj_dir, atlas_ref, "surf", "lh.white"
                )
                ctx_parc_lh_annot.append(tmp_annot)
                lh_obj = cltfree.AnnotParcellation.gii2annot(
                    gii_file=parc_file, annot_file=tmp_annot
                )

                tmp_annot = os.path.join(
                    fssubj_dir,
                    atlas_ref,
                    "label",
                    "rh." + at_name[0] + ".annot",
                )
                tmp_refsurf = os.path.join(
                    fssubj_dir, atlas_ref, "surf", "rh.white"
                )
                ctx_parc_rh_annot.append(tmp_annot)
                rh_obj = cltfree.AnnotParcellation.gii2annot(
                    gii_file=ctx_parc_rh[i], annot_file=tmp_annot
                )

        if not ctx_parc_lh_annot or not ctx_parc_rh_annot:
            raise ValueError(
                "Cortical parcellations should be supplied for both hemispheres."
            )
        else:
            meth_dict = {
                "method": method,
                "reference": atlas_ref,
                "labels": {
                    "lh": ctx_parc_lh_annot,
                    "rh": ctx_parc_rh_annot,
                },
            }

    elif atlas_src == "local":

        if atlas_type == "annot":
            atlas_dir = os.path.join(chim_dir, "data", "annot_atlases")
            atlas_ext = ".annot"
            method = "annot2indiv"

        elif atlas_type == "gcs":
            atlas_dir = os.path.join(chim_dir, "data", "gcs_atlases")
            atlas_ext = ".gcs"
            method = "gcs2indiv"
        else:
            raise ValueError(
                "The type of the cortical atlas is not valid. Please supply a valid type (annot or gcs)."
            )

        ctx_parc_lh = glob(os.path.join(atlas_dir, "*-L_*" + atlas_ext))
        ctx_parc_rh = glob(os.path.join(atlas_dir, "*-R_*" + atlas_ext))

        # Filtering for selecting the correct cortical parcellation
        ctx_parc_lh = cltmisc.filter_by_substring(
            ctx_parc_lh, atlas_names, bool_case=False
        )
        ctx_parc_rh = cltmisc.filter_by_substring(
            ctx_parc_rh, atlas_names, bool_case=False
        )
        ctx_parc_lh.sort()
        ctx_parc_rh.sort()

        if not ctx_parc_lh or not ctx_parc_rh:
            raise ValueError(
                "Cortical parcellations should be supplied for both hemispheres."
            )

        else:

            meth_dict = {
                "method": method,
                "reference": atlas_ref,
                "labels": {"lh": ctx_parc_lh, "rh": ctx_parc_rh},
            }

    return meth_dict


def _prepare_volumetric_templates(chim, supra, chim_dir, pipe_dict):
    """Prepare the volumetric templates for a (non-cortical) supra-region."""

    atlas_src = chim.parc_dict[supra]["source"]
    atlas_cad = chim.parc_dict[supra]["atlas"]
    type_cad = chim.parc_dict[supra]["type"]
    atlas_ref = chim.parc_dict[supra]["reference"]

    if atlas_src == "templateflow":

        # Getting the templates
        # Reference space
        tflow_home = _set_templateflow_home(
            pipe_dict["packages"]["templateflow"]["home_dir"]
        )
        t1_temp = tflow.get(
            atlas_ref,
            desc=None,
            resolution=[None, 1],
            suffix="T1w",
            extension="nii.gz",
        )

        # Getting the thalamic nuclei spams
        if type_cad == "spam":
            atlas_file = tflow.get(
                atlas_ref,
                desc=None,
                resolution=[None, 1],
                atlas=atlas_cad,
                suffix="probseg",
                extension="nii.gz",
            )

        elif type_cad == "maxprob":
            atlas_file = tflow.get(
                atlas_ref,
                desc=None,
                resolution=[None, 1],
                atlas=atlas_cad,
                suffix="dseg",
                extension="nii.gz",
            )
        else:
            # Raise an error if the type of the atlas is not valid and exit
            raise ValueError(
                "The type of the atlas is not valid. Please supply a valid type (spam or maxprob)."
            )

        meth_dict = {
            "method": "atlasbased",
            "type": type_cad,
            "reference": str(t1_temp),
            "labels": str(atlas_file),
        }

    elif atlas_src == "local":
        atlas_dir = os.path.join(chim_dir, "data", "vol_atlases")

        t1_temp = glob(os.path.join(atlas_dir, "*" + atlas_cad + "*_T1w.nii.gz"))

        if type_cad == "spam":
            atlas_file = glob(
                os.path.join(atlas_dir, "*" + atlas_cad + "*_probseg.nii.gz")
            )

        elif type_cad == "maxprob":
            atlas_file = glob(
                os.path.join(atlas_dir, "*" + atlas_cad + "*_dseg.nii.gz")
            )
        else:
            raise ValueError(
                "The type of the atlas is not valid. Please supply a valid type (spam or maxprob)."
            )

        meth_dict = {
            "method": "atlasbased",
            "type": type_cad,
            "reference": str(t1_temp),
            "labels": str(atlas_file),
        }

    elif atlas_src == "freesurfer":

        meth_dict = {
            "method": "comform2native",
            "type": None,
            "reference": "native",
            "labels": None,
        }

    elif atlas_src == "freesurferextra":
        meth_dict = {
            "method": "comform2native",
            "type": None,
            "reference": "native",
            "labels": atlas_src.lower(),
        }
    else:
        meth_dict = {
            "method": None,
            "type": None,
            "reference": "native",
            "labels": None,
        }

    return meth_dict
