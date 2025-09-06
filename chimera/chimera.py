#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from os import PathLike
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import bids
from bids import BIDSLayout
import nibabel as nib
from glob import glob
from operator import itemgetter
from datetime import datetime
import argparse
from typing import Union
import shutil
import uuid
import csv
import json
import subprocess
import scipy.ndimage as sc
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
from concurrent.futures import as_completed
import time
import copy
from threading import Lock

from templateflow import api as tflow


import clabtoolkit.misctools as cltmisc
from clabtoolkit.misctools import bcolors as bcolors
import clabtoolkit.freesurfertools as cltfree
import clabtoolkit.parcellationtools as cltparc
import clabtoolkit.bidstools as cltbids
import clabtoolkit.segmentationtools as cltseg
import clabtoolkit.imagetools as cltimg
from rich.progress import Progress

from .config_manager import (
    load_parcellations_info,
    _set_templateflow_home,
    _pipeline_info,
)

from .processing import launch_fsl_first

from .parcellation import (
    create_extra_regions_parc,
    _mix_side_prop,
    _print_availab_parcels,
)


# Define the Chimera class. This class will be used to create and work with Chimera objects
class Chimera:
    """
    Class to create and work with Chimera objects.

    Parameters
    ----------
    parc_code : str
        Parcellation code.

    parc_dict_file : str
        Parcellation dictionary file.

    supra_folder : str
        Folder containing the supraregions TSV files.

    Returns
    -------
    Chimera object


    """

    def __init__(
        self,
        parc_code: str,
        scale: Union[str, list] = None,
        seg: Union[str, list] = None,
        parc_dict_file: str = None,
        supra_folder: str = None,
    ):
        """
        Initialize the Chimera class

        Parameters
        ----------
        parc_code : str
            Parcellation code.

        parc_dict_file : str
            Parcellation dictionary file.
        supra_folder : str
            Folder containing the supraregions TSV files.

        Returns
        -------
        Chimera object

        Notes
        -----
        The parcellation code should be a string with the following format:
        "SFMIHISIFN", where each letter corresponds to a supra-region:

        """

        # Detecting the base directory
        chim_dir = os.path.dirname(os.path.abspath(__file__))

        # Rise an error if the parcellation code is not provided
        if parc_code is None:
            raise ValueError("Please provide a parcellation code")

        if parc_dict_file is not None:
            if not os.path.isfile(parc_dict_file):
                raise ValueError("The parcellation dictionary file does not exist")
        else:
            parc_dict_file = os.path.join(
                chim_dir, "config", "supraregions_dictionary.json"
            )

        if supra_folder is not None:
            if not os.path.isdir(supra_folder):
                raise ValueError(
                    "The the folder containing the supra-regions TSV files is not valid"
                )
            else:
                self.suprafolder = supra_folder
        else:
            self.suprafolder = os.path.join(chim_dir, "config", "supraregions")

        self.parc_dict, self.supra_dict = load_parcellations_info(
            parc_json=parc_dict_file, supra_folder=supra_folder
        )

        ####  Filtering the parcellation dictionary according to the parcellation code ####
        supra_names = list(self.parc_dict.keys())

        temp_dict = {}
        for i in range(len(parc_code)):
            if parc_code[i] in self.parc_dict[supra_names[i]].keys():

                # Defining the dictionary for the method
                meth_dict = {}
                meth_dict["code"] = parc_code[i]  # Add the parcellation code

                # Append the information of the parcellation using that method
                meth_dict.update(self.parc_dict[supra_names[i]][parc_code[i]])

                # Filtering the parcellation names by the scale and segmentation
                if i == 0:
                    parcel_names = meth_dict["parcels"]

                    # Filtering the parcellation names by the scale
                    if scale is not None:
                        if isinstance(scale, list):
                            # If the scale is a list do a loop over the elements and
                            # verify if the scale contains the string '_scale-'
                            scale_tmp = []
                            for sc in scale:
                                if "_scale-" not in sc:
                                    scale_tmp.append("_scale-" + sc)

                        elif isinstance(scale, str):
                            if "_scale-" not in scale:
                                scale_tmp = "_scale-" + scale

                        # Detect if the word scale is on any of the strings in parcel_names
                        if [s for s in parcel_names if "scale" in s]:
                            parcel_names = cltmisc.filter_by_substring(
                                parcel_names, scale_tmp, bool_case=False
                            )

                    # Filtering the parcellation names by the segmentation
                    if seg is not None:
                        if isinstance(seg, list):
                            # If the seg is a list do a loop over the elements and
                            # verify if the seg contains the string '_seg-'
                            seg_tmp = []
                            for sc in seg:
                                if "_seg-" not in sc:
                                    seg_tmp.append("_seg-" + sc)

                        elif isinstance(seg, str):
                            if "_seg-" not in seg:
                                seg_tmp = "_seg-" + seg

                        parcel_names = cltmisc.filter_by_substring(
                            parcel_names, seg_tmp, bool_case=False
                        )

                    # Saving the new parcels names
                    meth_dict["parcels"] = parcel_names

                # Adding the dictionary to the temp_dict
                temp_dict[supra_names[i]] = meth_dict

            else:
                # Print a message that the parcellation code is not present in the dictionary
                print(
                    "The parcellation code {} is not present in the dictionary for the supra-region {}.".format(
                        parc_code[i], supra_names[i]
                    )
                )
                # Error message and exit
                # sys.exit(1)

        self.parc_dict = temp_dict
        self.parc_code = parc_code
        self.scale = scale
        self.seg = seg

    def prepare_templates(self, fssubj_dir: str = None):
        """
        This method prepares the templates for the Chimera parcellation.
        Based on the code of the parcellation, it will download the necessary templates
        from the TemplateFlow repository or set up the templates directory for each supra region.
        It will also create the necessary symlinks to the FreeSurfer directory if the cortical
        parcellation is selected.

        Parameters
        ----------
        fssubj_dir : str
            FreeSurfer directory.

        Returns
        -------
        None

        Notes
        -----
        The FreeSurfer directory should be set up before running this method.
        The FreeSurfer directory should contain the subject used as reference for the cortical parcellation.


        Examples
        --------
        >>> chim = Chimera(parc_code="SFMIHISIFN", scale="100", seg="7n")
        >>> chim.prepare_templates(fssubj_dir="/path/to/freesurfer/subjects/fsaverage")


        """

        global pipe_dict

        # Setting up the FreeSurfer directory
        cltfree.FreeSurferSubject.set_freesurfer_directory(fssubj_dir)

        # Create the simlink to the FreeSurfer directory
        if "Cortical" in self.parc_dict.keys():
            cltfree.create_fsaverage_links(
                fssubj_dir,
                fsavg_dir=None,
                refsubj_name=self.parc_dict["Cortical"]["reference"],
            )

        # Detecting the base directory
        chim_dir = os.path.dirname(os.path.abspath(__file__))

        # Reading the names of the supra-regions
        supra_names = list(self.parc_dict.keys())

        for supra in supra_names:

            atlas_src = self.parc_dict[supra]["source"]
            atlas_str = self.parc_dict[supra]["atlas"]
            atlas_ref = self.parc_dict[supra]["reference"]

            if supra == "Cortical":

                # Atributes for the cortical parcellation
                atlas_type = self.parc_dict[supra]["type"]
                atlas_names = self.parc_dict[supra]["parcels"]

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
                                gii_file=parc_file,
                                ref_surf=tmp_refsurf,
                                annot_file=tmp_annot,
                                cont_tech=pipe_dict["packages"]["freesurfer"][
                                    "cont_tech"
                                ],
                                cont_image=pipe_dict["packages"]["freesurfer"][
                                    "container"
                                ],
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
                                gii_file=ctx_parc_rh[i],
                                ref_surf=tmp_refsurf,
                                annot_file=tmp_annot,
                                cont_tech=pipe_dict["packages"]["freesurfer"][
                                    "cont_tech"
                                ],
                                cont_image=pipe_dict["packages"]["freesurfer"][
                                    "container"
                                ],
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

            else:
                atlas_cad = self.parc_dict[supra]["atlas"]
                type_cad = self.parc_dict[supra]["type"]
                atlas_ref = self.parc_dict[supra]["reference"]

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

                    t1_temp = glob(
                        os.path.join(atlas_dir, "*" + atlas_cad + "*_T1w.nii.gz")
                    )

                    if type_cad == "spam":
                        atlas_file = glob(
                            os.path.join(
                                atlas_dir, "*" + atlas_cad + "*_probseg.nii.gz"
                            )
                        )

                    elif type_cad == "maxprob":
                        atlas_file = glob(
                            os.path.join(atlas_dir, "*" + atlas_cad + "*_dseg.nii.gz")
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

            self.parc_dict[supra]["processing"] = meth_dict

    def create_table(
        self,
        wm_index_offset: int = 3000,
        reg2rem: Union[list, str] = ["unknown", "medialwall", "corpuscallosum"],
    ):
        """
        This method creates a table with the regions that will be created using the Chimera parcellation.
        It allows to verify the region distribution for a specified parcelation code.

        Parameters
        ----------
        wm_index_offset : int
            Index offset for the white matter parcellation (default = 3000).

        reg2rem : list
            List of regions to remove from the parcellation.

        """

        # Detecting the base directory
        cwd = os.path.dirname(os.path.abspath(__file__))
        chim_dir = os.path.dirname(cwd)

        # Reading the names of the supra-regions
        supra_names = list(self.parc_dict.keys())

        lh_noctx_codes = []
        rh_noctx_codes = []
        lh_noctx_names = []
        rh_noctx_names = []
        lh_noctx_colors = []
        rh_noctx_colors = []
        bs_noctx_codes = []
        bs_noctx_names = []
        bs_noctx_colors = []

        ctx_parc_lh = []
        ctx_parc_rh = []

        desc_noctx = []
        for supra in supra_names:

            # Getting the information of the common atributes
            atlas_code = self.parc_dict[supra]["code"]
            atlas_str = self.parc_dict[supra]["atlas"]
            atlas_desc = self.parc_dict[supra]["description"]
            atlas_cita = self.parc_dict[supra]["citation"]
            atlas_src = self.parc_dict[supra]["source"]
            atlas_ref = self.parc_dict[supra]["reference"]

            if supra == "Cortical":

                # Atributes for the cortical parcellation
                atlas_type = self.parc_dict[supra]["type"]
                atlas_names = self.parc_dict[supra]["parcels"]

                # Selecting the source and downloading the parcellation
                if atlas_src == "templateflow":
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
                        ctx_parc_lh, or_filter=atlas_names, bool_case=False
                    )
                    ctx_parc_rh = cltmisc.filter_by_substring(
                        ctx_parc_rh, or_filter=atlas_names, bool_case=False
                    )

                elif atlas_src == "local":

                    if atlas_type == "annot":
                        atlas_dir = os.path.join(chim_dir, "data", "annot_atlases")
                        atlas_ext = ".annot"

                    elif atlas_type == "gcs":
                        atlas_dir = os.path.join(chim_dir, "data", "gcs_atlases")
                        atlas_ext = ".gcs"

                    ctx_parc_lh = glob(os.path.join(atlas_dir, "*-L_*" + atlas_ext))
                    ctx_parc_rh = glob(os.path.join(atlas_dir, "*-R_*" + atlas_ext))

                    # Filtering for selecting the correct cortical parcellation
                    ctx_parc_lh = cltmisc.filter_by_substring(
                        ctx_parc_lh, or_filter=atlas_names, bool_case=False
                    )
                    ctx_parc_rh = cltmisc.filter_by_substring(
                        ctx_parc_rh, or_filter=atlas_names, bool_case=False
                    )

            else:

                desc_noctx.append(atlas_desc)
                # Selecting the source and downloading the parcellation
                if atlas_src == "templateflow":

                    # Reference space
                    tflow_home = _set_templateflow_home(
                        pipe_dict["packages"]["templateflow"]["home_dir"]
                    )
                    ref_img = tflow.get(
                        atlas_ref,
                        desc=None,
                        resolution=1,
                        suffix="T1w",
                        extension="nii.gz",
                    )

                    # Getting the thalamic nuclei spams
                    parc_img = tflow.get(
                        atlas_ref,
                        desc=None,
                        resolution=1,
                        atlas=atlas_str,
                        suffix=atlas_type,
                        extension="nii.gz",
                    )

                if supra in self.supra_dict.keys():
                    meth_dict = self.parc_dict[supra]
                    st_dict = self.supra_dict[supra][supra][meth_dict["code"]]
                    if len(st_dict) == 1:
                        bs_noctx_codes = bs_noctx_codes + st_dict["mid"]["index"]
                        bs_noctx_names = bs_noctx_names + st_dict["mid"]["name"]
                        bs_noctx_colors = bs_noctx_colors + st_dict["mid"]["color"]

                    elif len(st_dict) == 2:
                        lh_noctx_codes = lh_noctx_codes + st_dict["lh"]["index"]
                        rh_noctx_codes = rh_noctx_codes + st_dict["rh"]["index"]

                        lh_noctx_names = lh_noctx_names + st_dict["lh"]["name"]
                        rh_noctx_names = rh_noctx_names + st_dict["rh"]["name"]

                        lh_noctx_colors = lh_noctx_colors + st_dict["lh"]["color"]
                        rh_noctx_colors = rh_noctx_colors + st_dict["rh"]["color"]

        if rh_noctx_names:
            indexes = cltmisc.get_indexes_by_substring(rh_noctx_names, reg2rem).tolist()
            # Remove the elements in all_names and all_colors
            if indexes:
                for i in indexes:
                    rh_noctx_names.pop(i)
                    rh_noctx_codes.pop(i)
                    rh_noctx_colors.pop(i)

        if lh_noctx_names:
            indexes = cltmisc.get_indexes_by_substring(lh_noctx_names, reg2rem).tolist()
            # Remove the elements in all_names and all_colors
            if indexes:
                for i in indexes:
                    lh_noctx_names.pop(i)
                    lh_noctx_codes.pop(i)
                    lh_noctx_colors.pop(i)

        if bs_noctx_names:
            indexes = cltmisc.get_indexes_by_substring(bs_noctx_names, reg2rem).tolist()
            # Remove the elements in all_names and all_colors
            if indexes:
                for i in indexes:
                    bs_noctx_names.pop(i)
                    bs_noctx_codes.pop(i)
                    bs_noctx_codes.pop(i)

        # Creating the list of dataframes for the different parcellations
        tab_list = []
        desc_list = []
        parc_id = "atlas-chimera" + self.parc_code
        parc_id_list = []

        # If ctx_parc_lh is empty, it means that the parcellation is not available
        if len(ctx_parc_lh) == 0:
            all_names = rh_noctx_names + lh_noctx_names + bs_noctx_names
            all_colors = rh_noctx_colors + lh_noctx_colors + bs_noctx_colors
            index = np.arange(1, len(all_names) + 1).tolist()
            tab_df = pd.DataFrame(
                {"index": index, "name": all_names, "color": all_colors}
            )
            tab_list.append(tab_df)

            gen_desc = ["# Parcellation code: " + self.parc_code]
            gen_desc.append(desc_noctx)

            parc_id_list.append(parc_id)

        else:
            for i, parc_file in enumerate(ctx_parc_lh):

                gen_desc = ["# Parcellation code: " + self.parc_code]

                tmp_name = os.path.basename(ctx_parc_lh[i])
                tmp_ent = tmp_name.split("_")[:-1]

                # Get the element that contains the string 'scale' and extract it
                scale_ent = [s for s in tmp_ent if "scale" in s]
                if scale_ent:
                    scale_ent = scale_ent[0]
                    scale_ent = scale_ent.split("-")[1]
                    parc_id = parc_id + "_scale-" + scale_ent

                    # Add a the segmentation to to the string of the general description list
                    gen_desc[0] = gen_desc[0] + ". Scale: " + scale_ent

                # Get the element that contains the string 'seg' and extract it
                seg_ent = [s for s in tmp_ent if "seg" in s]

                if seg_ent:
                    seg_ent = seg_ent[0]
                    seg_ent = seg_ent.split("-")[1]
                    parc_id = parc_id + "_seg-" + seg_ent

                    # Add a the segmentation to to the string of the general description list
                    gen_desc[0] = gen_desc[0] + ". Segmentation: " + seg_ent

                gen_desc.append(self.parc_dict["Cortical"]["description"])
                gen_desc = gen_desc + desc_noctx

                # Reading the cortical parcellations
                lh_obj = cltfree.AnnotParcellation(parc_file=ctx_parc_lh[i])
                rh_obj = cltfree.AnnotParcellation(parc_file=ctx_parc_rh[i])

                df_lh, out_tsv = lh_obj.export_to_tsv(prefix2add="ctx-lh-")
                df_rh, out_tsv = rh_obj.export_to_tsv(prefix2add="ctx-rh-")

                # Convert the column name of the dataframe to a list
                lh_ctx_name = df_lh["name"].tolist()
                rh_ctx_name = df_rh["name"].tolist()

                # Convert the column color of the dataframe to a list
                lh_ctx_color = df_lh["color"].tolist()
                rh_ctx_color = df_rh["color"].tolist()

                ## Removing elements from the table according to their name for both
                indexes = cltmisc.get_indexes_by_substring(
                    lh_ctx_name, reg2rem
                ).tolist()
                if indexes:
                    for i in indexes:
                        lh_ctx_name.pop(i)
                        lh_ctx_color.pop(i)

                indexes = cltmisc.get_indexes_by_substring(
                    rh_ctx_name, reg2rem
                ).tolist()
                if indexes:
                    for i in indexes:
                        rh_ctx_name.pop(i)
                        rh_ctx_color.pop(i)

                # Concatenating the lists
                if "GyralWM" in self.parc_dict.keys():
                    gen_desc.append(self.parc_dict["GyralWM"]["description"])

                    wm_rh_name = cltmisc.correct_names(
                        rh_ctx_name, replace=["ctx-rh-", "wm-rh-"]
                    )
                    wm_rh_indexes = np.arange(1, len(wm_rh_name) + 1) + wm_index_offset
                    wm_rh_indexes = wm_rh_indexes.tolist()

                    wm_lh_name = cltmisc.correct_names(
                        lh_ctx_name, replace=["ctx-lh-", "wm-lh-"]
                    )
                    wm_lh_indexes = (
                        np.arange(1, len(wm_lh_name) + 1)
                        + len(rh_ctx_name)
                        + len(rh_noctx_names)
                        + wm_index_offset
                    )
                    wm_lh_indexes = wm_lh_indexes.tolist()

                    wm_rh_color = rh_ctx_color
                    wm_lh_color = lh_ctx_color

                else:
                    wm_lh_indexes = []
                    wm_rh_indexes = []
                    wm_lh_name = []
                    wm_rh_name = []
                    wm_lh_color = []
                    wm_rh_color = []

                # Right hemisphere
                rh_all_names = rh_ctx_name + rh_noctx_names
                rh_all_indexes = np.arange(1, len(rh_all_names) + 1).tolist()

                # Left hemisphere
                lh_all_names = (
                    lh_ctx_name
                    + lh_noctx_names
                    + bs_noctx_names
                    + wm_rh_name
                    + wm_lh_name
                )
                lh_all_indexes = np.arange(
                    1, len(lh_ctx_name + lh_noctx_names + bs_noctx_names) + 1
                ) + np.max(rh_all_indexes)
                lh_all_indexes = lh_all_indexes.tolist()

                rh_all_colors = rh_ctx_color + rh_noctx_colors
                lh_all_colors = (
                    lh_ctx_color
                    + lh_noctx_colors
                    + bs_noctx_colors
                    + wm_rh_color
                    + wm_lh_color
                )

                # Concatenating the hemispheres
                all_names = rh_all_names + lh_all_names
                all_colors = rh_all_colors + lh_all_colors
                all_indexes = (
                    rh_all_indexes + lh_all_indexes + wm_rh_indexes + wm_lh_indexes
                )

                # Generating a dataframe
                tab_df = pd.DataFrame(
                    {"index": all_indexes, "name": all_names, "color": all_colors}
                )
                tab_list.append(tab_df)
                desc_list.append(gen_desc)
                parc_id_list.append(parc_id)

        # Add the tab_list as an attribute of the class
        self.regtable = {"parc_id": parc_id_list, "desc": desc_list, "table": tab_list}

    def export_table(self, out_basename: str = None, format: Union[list, str] = "tsv"):
        """
        This method exports the table of the regions to a TSV or a LUT file.

        Parameters
        ----------
        out_basename : str
            Output basename for the TSV or LUT file.

        format : str or list
            Format of the output file. It can be 'tsv', 'lut' or ['tsv, lut'].

        """

        if out_basename is None:
            # Exit if the output basename is not provided
            raise ValueError(
                "Please provide an output basename for the TSV or LUT file."
            )

        out_name = os.path.basename(out_basename)

        out_dir = os.path.dirname(out_basename)
        out_dir = Path(out_dir)

        # Create the output directory if it does not exist
        out_dir.mkdir(parents=True, exist_ok=True)

        parc_ids = self.regtable["parc_id"]
        parc_desc = self.regtable["desc"]
        parc_tables = self.regtable["table"]

        # Export the table to a TSV file
        for i, tab_df in enumerate(parc_tables):

            if isinstance(format, list):
                if "tsv" in format:
                    out_file_tsv = os.path.join(
                        str(out_dir), out_name + "_" + parc_ids[i] + ".tsv"
                    )
                    cltparc.Parcellation.write_tsvtable(
                        tsv_df=tab_df, out_file=out_file_tsv, force=True
                    )

                if "lut" in format:
                    out_file_lut = os.path.join(
                        str(out_dir), out_name + "_" + parc_ids[i] + ".lut"
                    )
                    codes = tab_df["index"].tolist()
                    names = tab_df["name"].tolist()
                    colors = tab_df["color"].tolist()
                    cltparc.Parcellation.write_luttable(
                        codes=codes,
                        names=names,
                        colors=colors,
                        out_file=out_file_lut,
                        headerlines=parc_desc,
                        force=True,
                    )
            else:
                if format == "tsv":
                    out_file_tsv = os.path.join(
                        str(out_dir), out_name + "_" + parc_ids[i] + ".tsv"
                    )
                    cltparc.Parcellation.write_tsvtable(
                        tsv_df=tab_df, out_file=out_file_tsv, force=True
                    )

                if format == "lut":
                    out_file_lut = os.path.join(
                        str(out_dir), out_name + "_" + parc_ids[i] + ".lut"
                    )
                    codes = tab_df["index"].tolist()
                    names = tab_df["name"].tolist()
                    colors = tab_df["color"].tolist()
                    cltparc.Parcellation.write_luttable(
                        codes=codes,
                        names=names,
                        colors=colors,
                        out_file=out_file_lut,
                        headerlines=parc_desc,
                        force=True,
                    )

    def build_lut_header(self):
        """
        This method builds the header of the LUT file.

        """

        # Detecting the base directory
        chim_dir = os.path.dirname(os.path.abspath(__file__))

        # Get the absolute of this file
        parc_json = os.path.join(chim_dir, "config", "supraregions_dictionary.json")

        # Reading the parcellation dictionary
        with open(parc_json) as f:
            parc_dict = json.load(f)

        # Reading the names of the supra-regions
        supra_names = list(parc_dict.keys())

        # Reading the parcellation code
        chim_code = self.parc_code

        # Creating the header lines
        headerlines = [" # Chimera parcellation code: {}".format(self.parc_code)]
        for i, supra in enumerate(supra_names):
            tmp_dict = parc_dict[supra]

            # Check if the parcellation code is in the dictionary
            if chim_code[i] in tmp_dict.keys():
                tmp_dict = tmp_dict[chim_code[i]]

                if tmp_dict["description"].endswith("."):
                    tmp_dict["description"] = tmp_dict["description"][:-1]

                cite = "{} {}.".format(tmp_dict["atlas"], tmp_dict["citation"])
                glob_desc = tmp_dict["description"] + ". Name: " + cite
                headerlines.append("    " + glob_desc)
            else:
                headerlines.append(
                    "    # {}. The parcellation code {} is not present in the dictionary for the supra-region {}.".format(
                        i + 1, chim_code[i], supra
                    )
                )

        return headerlines

    def build_parcellation(
        self,
        t1: str,
        bids_dir: str,
        deriv_dir: str = None,
        fssubj_dir: str = None,
        growwm: Union[str, int] = None,
        bool_mixwm: bool = False,
        force: bool = False,
    ):
        """
        This method builds the parcellation for the selected parcellation code.

        Parameters
        ----------
        t1 : str
            T1-weighted image in the BIDs format.

        bids_dir : str
            BIDs dataset directory.

        deriv_dir : str
            BIDs derivative directory.

        fssubj_dir : str
            FreeSurfer subjects directory. If not provided, the environment variable FREESURFER_HOME
            will be used.

        growwm : str or int
            Grow of GM labels inside the white matter in mm.
            If None, no growing will be applied.
            If "0", no growing will be applied.
            If int, the value of growing in mm.
            If a list of int, the values of growing in mm.

        bool_mixwm : bool
            If True, the growing will be applied to the voxels that are in contact with the white matter
            and they will mix with the cortical labels.

        force : bool
            Overwrite the results.

        """
        global pipe_dict

        # Detecting the base directory
        cwd = os.path.dirname(os.path.abspath(__file__))
        chim_dir = os.path.dirname(cwd)

        if not os.path.isfile(t1):
            raise ValueError("Please provide a valid T1 image.")

        # Getting the entities from the name
        anat_dir = os.path.dirname(t1)
        t1_name = os.path.basename(t1)
        ent_dict = cltbids.str2entity(t1_name)

        temp_entities = t1_name.split("_")[:-1]
        fullid = "_".join(temp_entities)
        ent_dict_fullid = cltbids.str2entity(fullid)

        if "ses" in ent_dict.keys():
            path_cad = "sub-" + ent_dict["sub"] + os.path.sep + "ses-" + ent_dict["ses"]
        else:
            path_cad = "sub-" + ent_dict["sub"]

        # Creating Chimera directories
        if deriv_dir is None:
            chim_dir = os.path.join(
                bids_dir, "derivatives", "chimera", path_cad, "anat"
            )
        else:
            chim_dir = os.path.join(deriv_dir, "chimera", path_cad, "anat")

        # Create the Chimera directory if it does not exist
        chim_dir = Path(chim_dir)
        chim_dir.mkdir(parents=True, exist_ok=True)

        supra_names = list(self.parc_dict.keys())
        # Detecting if Cortical is on the list of supra-regions
        bool_ctx = False
        if "Cortical" in supra_names:
            # Remove it from the list
            supra_names.remove("Cortical")
            bool_ctx = True

        # ----------- Veryfing the existence of the parcellations, otherwise, compute them  --------- #
        if force:
            bool_chim_exist = False
        else:
            bool_chim_exist = True

            if bool_ctx:

                # Atributes for the cortical parcellation
                atlas_names = self.parc_dict["Cortical"]["parcels"]
                for at_name in atlas_names:
                    ## -------- Cortical parcellation for the left hemisphere ---------------
                    # Creating the name for the output file

                    if growwm is None:
                        growwm = ["0"]

                    for ngrow in np.arange(len(growwm)):
                        if growwm[ngrow] == "0":
                            out_vol_name = fullid + "_" + at_name + "_dseg.nii.gz"
                        else:

                            ent_dict = cltbids.str2entity(at_name)
                            if "desc" in ent_dict.keys():
                                if bool_mixwm:
                                    ent_dict["desc"] = (
                                        ent_dict["desc"]
                                        + "grow"
                                        + str(growwm[ngrow])
                                        + "mm+mixwm"
                                    )
                                else:
                                    ent_dict["desc"] = (
                                        ent_dict["desc"]
                                        + "grow"
                                        + str(growwm[ngrow])
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
                                        + str(growwm[ngrow])
                                        + "mm+mixwm_dseg.nii.gz"
                                    )
                                else:
                                    out_vol_name = (
                                        fullid
                                        + "_"
                                        + at_name
                                        + "_desc-grow"
                                        + str(growwm[ngrow])
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

        # ------- End of veryfing the existence of the parcellations, otherwise, compute them  --------- #

        # Run chimera if the desired parcellation does not exist
        if not bool_chim_exist:

            ######## ----------- Detecting FREESURFER_HOME directory ------------- #
            if pipe_dict["packages"]["freesurfer"]["cont_tech"] != "local":
                cont_tech = pipe_dict["packages"]["freesurfer"]["cont_tech"]
                cont_image = pipe_dict["packages"]["freesurfer"]["container"]

                # Detecting the FreeSurfer home directory using the container
                cmd_bashargs = ["echo", "$FREESURFER_HOME"]
                cmd_cont = cltmisc.generate_container_command(
                    cmd_bashargs, cont_tech, cont_image
                )  # Generating container command
                out_cmd = subprocess.run(
                    cmd_cont, stdout=subprocess.PIPE, universal_newlines=True
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
                subprocess.run(
                    cmd_cont, stdout=subprocess.PIPE, universal_newlines=True
                )
                fslut_file = os.path.join("/tmp", tmp_name)

                ######## ------------- Reading FreeSurfer color lut table ------------ #
                lut_dict = cltparc.Parcellation.read_luttable(fslut_file)

                os.remove(fslut_file)

            else:

                fshome_dir = os.getenv("FREESURFER_HOME")
                fslut_file = os.path.join(fshome_dir, "FreeSurferColorLUT.txt")

                ######## ------------- Reading FreeSurfer color lut table ------------ #
                lut_dict = cltparc.Parcellation.read_luttable(fslut_file)

            # Extracting the information from the lut dictionary
            st_codes = lut_dict["index"]
            st_names = lut_dict["name"]
            st_colors = lut_dict["color"]

            ######## ----- Running FreeSurfer if it was not previously computed ------ #
            sub2proc = cltfree.FreeSurferSubject(fullid, subjs_dir=fssubj_dir)

            cont_tech_freesurfer = pipe_dict["packages"]["freesurfer"]["cont_tech"]
            cont_image_freesurfer = pipe_dict["packages"]["freesurfer"]["container"]
            freesurfer_license_file = pipe_dict["packages"]["freesurfer"]["license"]
            cont_tech_ants = pipe_dict["packages"]["ants"]["cont_tech"]
            cont_image_ants = pipe_dict["packages"]["ants"]["container"]
            cont_tech_fsl = pipe_dict["packages"]["fsl"]["cont_tech"]
            cont_image_fsl = pipe_dict["packages"]["fsl"]["container"]

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

            if "aseg_parc" not in locals():
                aseg_parc = cltparc.Parcellation(parc_file=nii_image)
                aseg_parc.index = st_codes
                aseg_parc.name = st_names
                aseg_parc.color = st_colors
                aseg_parc.adjust_values()

            # Creating the parcellation for the extra regions
            extra_parc = create_extra_regions_parc(aparc=nii_image)

            # Remove the nifti file
            os.remove(nii_image)

            # Building the main header information for the LUT file
            glob_header_info = self.build_lut_header()

            # Processing that will be perfomed for multiple supra-regions
            gm_sub_names = list(self.parc_dict.keys())

            # Remove 'Cortical', 'GyralWM' and 'WhiteMatter' from the list
            if "Cortical" in gm_sub_names:
                gm_sub_names.remove("Cortical")

            if "GyralWM" in gm_sub_names:
                gm_sub_names.remove("GyralWM")

            if "WhiteMatter" in supra_names:
                gm_sub_names.remove("WhiteMatter")

            # Taking the image dimensions and the affine matrix in native space
            t1_image = nib.load(t1)
            affine = t1_image.affine

            # Create a numpy array with the same dimensions of the T1 image and fill it with zeros.
            # The array will be used to store the parcellation. The elements should be integers.
            ref_image = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
            lh2refill = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
            rh2refill = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)

            # Creating the parcellation objects
            lh_parc = cltparc.Parcellation(
                parc_file=ref_image, affine=affine
            )  # It will include the parcellation for the left hemisphere
            rh_parc = copy.deepcopy(
                lh_parc
            )  # It will include the parcellation for the right hemisphere
            mid_parc = copy.deepcopy(
                lh_parc
            )  # It will include the parcellation for structures that do not belong to any hemisphere

            files2del = []  # Temporal files that will be deleted
            exec_cmds = []
            for supra in gm_sub_names:

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
                if atlas_type == "spam":
                    # Verifying the existence of the threshold
                    if "probthresh" not in self.parc_dict[supra].keys():
                        spam_thresh = 0.05
                    else:
                        spam_thresh = self.parc_dict[supra]["probthresh"]

                deriv_fold = self.parc_dict[supra]["deriv_volfold"]
                proc_dict = self.parc_dict[supra]["processing"]

                if proc_dict["method"] == "comform2native":

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
                            raise ValueError(
                                "The Freesurfer extra parcellation was not found."
                            )

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
                            lh_supra_parc.index = self.supra_dict[supra][supra][
                                atlas_code
                            ]["lh"]["index"]
                            lh_supra_parc.name = self.supra_dict[supra][supra][
                                atlas_code
                            ]["lh"]["name"]
                            lh_supra_parc.color = self.supra_dict[supra][supra][
                                atlas_code
                            ]["lh"]["color"]
                            lh_supra_parc.keep_by_code(
                                codes2keep=self.supra_dict[supra][supra][atlas_code][
                                    "lh"
                                ]["index"]
                            )
                            lh_supra_parc.export_colortable(
                                out_file=lh_nii_image.replace(".nii.gz", ".lut"),
                                lut_type="lut",
                            )
                            lh_supra_parc.export_colortable(
                                out_file=lh_nii_image.replace(".nii.gz", ".tsv"),
                                lut_type="tsv",
                            )

                            rh_supra_parc = cltparc.Parcellation(parc_file=rh_nii_image)
                            rh_supra_parc.index = self.supra_dict[supra][supra][
                                atlas_code
                            ]["rh"]["index"]
                            rh_supra_parc.name = self.supra_dict[supra][supra][
                                atlas_code
                            ]["rh"]["name"]
                            rh_supra_parc.color = self.supra_dict[supra][supra][
                                atlas_code
                            ]["rh"]["color"]
                            rh_supra_parc.keep_by_code(
                                codes2keep=self.supra_dict[supra][supra][atlas_code][
                                    "rh"
                                ]["index"]
                            )
                            rh_supra_parc.export_colortable(
                                out_file=rh_nii_image.replace(".nii.gz", ".lut"),
                                lut_type="lut",
                            )
                            rh_supra_parc.export_colortable(
                                out_file=rh_nii_image.replace(".nii.gz", ".tsv"),
                                lut_type="tsv",
                            )

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
                            index, name, color = _mix_side_prop(
                                self.supra_dict[supra][supra][atlas_code]
                            )
                            tmp_parc.index = index
                            tmp_parc.name = name
                            tmp_parc.color = color

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
                                lh_supra_parc = copy.deepcopy(tmp_parc)
                                lh_supra_parc.keep_by_code(
                                    codes2keep=self.supra_dict[supra][supra][
                                        atlas_code
                                    ]["lh"]["index"]
                                )

                            # Right Hemisphere
                            if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
                                rh_supra_parc = copy.deepcopy(tmp_parc)
                                rh_supra_parc.keep_by_code(
                                    codes2keep=self.supra_dict[supra][supra][
                                        atlas_code
                                    ]["rh"]["index"]
                                )

                            # Non-hemispheric structures
                            if (
                                "mid"
                                in self.supra_dict[supra][supra][atlas_code].keys()
                            ):
                                mid_supra_parc = copy.deepcopy(tmp_parc)
                                mid_supra_parc.keep_by_code(
                                    codes2keep=self.supra_dict[supra][supra][
                                        atlas_code
                                    ]["mid"]["index"]
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

                        if "aseg_parc" not in locals():
                            if not os.path.isfile(nii_image):
                                sub2proc.conform2native(
                                    mgz_conform=mgz_image,
                                    nii_native=nii_image,
                                    cont_image=cont_image_freesurfer,
                                    cont_tech=cont_tech_freesurfer,
                                    force=force,
                                )
                                files2del.append(nii_image)

                            tmp_parc = cltparc.Parcellation(parc_file=nii_image)

                        else:
                            tmp_parc = copy.deepcopy(aseg_parc)

                        # Left Hemisphere
                        if "lh" in self.supra_dict[supra][supra][atlas_code].keys():
                            lh_supra_parc = copy.deepcopy(tmp_parc)
                            lh_supra_parc.index = self.supra_dict[supra][supra][
                                atlas_code
                            ]["lh"]["index"]
                            lh_supra_parc.name = self.supra_dict[supra][supra][
                                atlas_code
                            ]["lh"]["name"]
                            lh_supra_parc.color = self.supra_dict[supra][supra][
                                atlas_code
                            ]["lh"]["color"]
                            lh_supra_parc.keep_by_code(
                                codes2keep=self.supra_dict[supra][supra][atlas_code][
                                    "lh"
                                ]["index"]
                            )

                        # Right Hemisphere
                        if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
                            rh_supra_parc = copy.deepcopy(tmp_parc)
                            rh_supra_parc.index = self.supra_dict[supra][supra][
                                atlas_code
                            ]["rh"]["index"]
                            rh_supra_parc.name = self.supra_dict[supra][supra][
                                atlas_code
                            ]["rh"]["name"]
                            rh_supra_parc.color = self.supra_dict[supra][supra][
                                atlas_code
                            ]["rh"]["color"]
                            rh_supra_parc.keep_by_code(
                                codes2keep=self.supra_dict[supra][supra][atlas_code][
                                    "rh"
                                ]["index"]
                            )

                        # Non-hemispheric structures
                        if "mid" in self.supra_dict[supra][supra][atlas_code].keys():
                            mid_supra_parc = copy.deepcopy(tmp_parc)
                            mid_supra_parc.index = self.supra_dict[supra][supra][
                                atlas_code
                            ]["mid"]["index"]
                            mid_supra_parc.name = self.supra_dict[supra][supra][
                                atlas_code
                            ]["mid"]["name"]
                            mid_supra_parc.color = self.supra_dict[supra][supra][
                                atlas_code
                            ]["mid"]["color"]
                            mid_supra_parc.keep_by_code(
                                codes2keep=self.supra_dict[supra][supra][atlas_code][
                                    "mid"
                                ]["index"]
                            )

                elif proc_dict["method"] == None:

                    # Running FIRST if it is needed
                    if atlas_code == "R":
                        fsl_outdir = os.path.join(
                            deriv_dir, deriv_fold, path_cad, "anat"
                        )
                        first_nii = os.path.join(
                            str(fsl_outdir),
                            fullid + "_atlas-" + atlas_str + "_dseg.nii.gz",
                        )

                        if "first_parc" not in locals():
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
                            first_parc = cltparc.Parcellation(parc_file=first_nii)

                        tmp_parc = copy.deepcopy(first_parc)
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
                        lh_supra_parc = copy.deepcopy(tmp_parc)
                        lh_supra_parc.keep_by_code(
                            codes2keep=self.supra_dict[supra][supra][atlas_code]["lh"][
                                "index"
                            ]
                        )

                    # Right Hemisphere
                    if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
                        rh_supra_parc = copy.deepcopy(tmp_parc)
                        rh_supra_parc.keep_by_code(
                            codes2keep=self.supra_dict[supra][supra][atlas_code]["rh"][
                                "index"
                            ]
                        )

                    # Non-hemispheric structures
                    if "mid" in self.supra_dict[supra][supra][atlas_code].keys():
                        mid_supra_parc = copy.deepcopy(tmp_parc)
                        mid_supra_parc.keep_by_code(
                            codes2keep=self.supra_dict[supra][supra][atlas_code]["mid"][
                                "index"
                            ]
                        )

                elif proc_dict["method"] == "atlasbased":
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
                        str(work_dir),
                        fullid + "_atlas-" + atlas_str + "_probseg.nii.gz",
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
                        or force
                    ):
                        work_dir.mkdir(parents=True, exist_ok=True)

                        # Detecting the side
                        sides_ids = list(
                            self.supra_dict[supra][supra][atlas_code].keys()
                        )
                        sides_ids = sorted(
                            sides_ids, key=lambda x: not ("lh" in x or "rh" in x)
                        )

                        # Masking the cerebellum from T1w image
                        tmp_t1 = t1
                        if supra == "Cerebellum":
                            if self.parc_dict[supra]["name"] == "SUIT":
                                tmp_t1 = os.path.join(
                                    str(work_dir), "tmp_cerb_bs.nii.gz"
                                )
                                files2del.append(tmp_t1)

                                cltimg.crop_image_from_mask(
                                    t1,
                                    aseg_parc.data,
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
                                    cltimg.cropped_to_native(
                                        out_parc_spam, t1, out_parc_spam
                                    )

                            for side_cont, side in enumerate(sides_ids):
                                vol_indexes = (
                                    np.array(
                                        self.supra_dict[supra][supra][atlas_code][side][
                                            "index"
                                        ]
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
                                files2del.append(tmp_par_file)

                                tmp_parc_file = cltimg.spams2maxprob(
                                    out_parc_spam,
                                    prob_thresh=spam_thresh,
                                    vol_indexes=vol_indexes,
                                    maxp_name=tmp_par_file,
                                )

                                tmp_parc = cltparc.Parcellation(parc_file=tmp_parc_file)
                                tmp_parc.index = vol_indexes + 1
                                tmp_parc.name = self.supra_dict[supra][supra][
                                    atlas_code
                                ][side]["name"]
                                tmp_parc.color = self.supra_dict[supra][supra][
                                    atlas_code
                                ][side]["color"]

                                if side in self.supra_dict[supra][supra]["F"].keys():
                                    aseg_code = self.supra_dict[supra][supra]["F"][
                                        side
                                    ]["index"]
                                    tmp_parc.apply_mask(
                                        image_mask=aseg_parc,
                                        codes2mask=aseg_code,
                                        fill=True,
                                    )

                                else:
                                    # Selecting all the region in case there is no definition of left and right hemispheres
                                    all_reg_codes = []
                                    for side_f in self.supra_dict[supra][supra][
                                        "F"
                                    ].keys():
                                        aseg_code = self.supra_dict[supra][supra]["F"][
                                            side_f
                                        ]["index"]
                                        all_reg_codes = all_reg_codes + aseg_code

                                    glob_mask_parc = copy.deepcopy(aseg_parc)
                                    glob_mask_parc.group_by_code(
                                        codes2group=all_reg_codes, new_codes=1
                                    )
                                    tmp_parc.apply_mask(
                                        image_mask=glob_mask_parc, codes2mask=1
                                    )

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
                                save_lut=True,
                                save_tsv=True,
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
                                    cltimg.cropped_to_native(
                                        out_parc_maxp, t1, out_parc_maxp
                                    )

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
                                files2del.append(tmp_par_file)

                                tmp_parc = cltparc.Parcellation(parc_file=out_parc_maxp)
                                tmp_parc.data = np.round(tmp_parc.data)
                                tmp_parc.index = np.array(
                                    self.supra_dict[supra][supra][atlas_code][side][
                                        "index"
                                    ]
                                )
                                tmp_parc.name = self.supra_dict[supra][supra][
                                    atlas_code
                                ][side]["name"]
                                tmp_parc.color = self.supra_dict[supra][supra][
                                    atlas_code
                                ][side]["color"]

                                if side in self.supra_dict[supra][supra]["F"].keys():
                                    aseg_code = self.supra_dict[supra][supra]["F"][
                                        side
                                    ]["index"]
                                    tmp_parc.apply_mask(
                                        image_mask=aseg_parc,
                                        codes2mask=aseg_code,
                                        fill=True,
                                    )

                                else:
                                    # Selecting all the region in case there is no definition of left and right hemispheres
                                    all_reg_codes = []
                                    for side_f in self.supra_dict[supra][supra][
                                        "F"
                                    ].keys():
                                        aseg_code = self.supra_dict[supra][supra]["F"][
                                            side_f
                                        ]["index"]
                                        all_reg_codes = all_reg_codes + aseg_code

                                    glob_mask_parc = copy.deepcopy(aseg_parc)
                                    glob_mask_parc.group_by_code(
                                        codes2group=all_reg_codes, new_codes=1
                                    )
                                    tmp_parc.apply_mask(
                                        image_mask=glob_mask_parc, codes2mask=1
                                    )

                                # Adjusting the values to the ones existing on the 3D image
                                tmp_parc.adjust_values()
                                if side_cont == 0:
                                    def_parc = copy.deepcopy(tmp_parc)
                                else:
                                    def_parc.add_parcellation(tmp_parc)

                            def_parc.save_parcellation(
                                out_file=out_parc_maxp,
                                affine=def_parc.affine,
                                save_lut=True,
                                save_tsv=True,
                            )

                    tmp_parc = cltparc.Parcellation(parc_file=out_parc_maxp)
                    index, name, color = _mix_side_prop(
                        self.supra_dict[supra][supra][atlas_code]
                    )
                    tmp_parc.index = index
                    tmp_parc.name = name
                    tmp_parc.color = color
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
                        lh_supra_parc = copy.deepcopy(tmp_parc)
                        lh_supra_parc.keep_by_code(
                            codes2keep=self.supra_dict[supra][supra][atlas_code]["lh"][
                                "index"
                            ]
                        )

                    # Right Hemisphere
                    if "rh" in self.supra_dict[supra][supra][atlas_code].keys():
                        rh_supra_parc = copy.deepcopy(tmp_parc)
                        rh_supra_parc.keep_by_code(
                            codes2keep=self.supra_dict[supra][supra][atlas_code]["rh"][
                                "index"
                            ]
                        )

                    # Non-hemispheric structures
                    if "mid" in self.supra_dict[supra][supra][atlas_code].keys():
                        mid_supra_parc = copy.deepcopy(tmp_parc)
                        mid_supra_parc.keep_by_code(
                            codes2keep=self.supra_dict[supra][supra][atlas_code]["mid"][
                                "index"
                            ]
                        )

                # Appending the parcellations
                if "lh_supra_parc" in locals():
                    lh_supra_parc.rearrange_parc()

                    if "F" in self.supra_dict[supra][supra].keys():
                        # Use the FreeSurfer parcellation to detect the voxels that are not in the lh_supra_parc
                        lh2refill_parc = copy.deepcopy(aseg_parc)
                        lh2refill_parc.keep_by_code(
                            codes2keep=self.supra_dict[supra][supra]["F"]["lh"]["index"]
                        )

                        # Find the voxels that are not in the lh_supra_parc and are in the lh2refill
                        ind = np.where(
                            (lh_supra_parc.data == 0) & (lh2refill_parc.data != 0)
                        )
                        lh2refill[ind] = 1
                        del lh2refill_parc

                    # Add the parcellation to the global left subcortical parcellation
                    lh_parc.add_parcellation(lh_supra_parc, append=True)
                    nlh_subc = len(lh_parc.index)
                    del lh_supra_parc
                    # lh_parc.save_parcellation(out_file= '/home/yaleman/lh_test.nii.gz', save_lut=True)

                if "rh_supra_parc" in locals():
                    rh_supra_parc.rearrange_parc()

                    if "F" in self.supra_dict[supra][supra].keys():
                        # Use the FreeSurfer parcellation to detect the voxels that are not in the lh_supra_parc
                        rh2refill_parc = copy.deepcopy(aseg_parc)
                        rh2refill_parc.keep_by_code(
                            codes2keep=self.supra_dict[supra][supra]["F"]["rh"]["index"]
                        )

                        # Find the voxels that are not in the lh_supra_parc and are in the lh2refill
                        ind = np.where(
                            (rh_supra_parc.data == 0) & (rh2refill_parc.data != 0)
                        )
                        rh2refill[ind] = 1
                        del rh2refill_parc

                    # Add the parcellation to the global right subcortical parcellation
                    rh_parc.add_parcellation(rh_supra_parc, append=True)
                    nrh_subc = len(rh_parc.index)
                    del rh_supra_parc
                    # rh_parc.save_parcellation(out_file= '/home/yaleman/rh_test.nii.gz', save_lut=True)

                if "mid_supra_parc" in locals():
                    mid_supra_parc.rearrange_parc()
                    mid_parc.add_parcellation(mid_supra_parc, append=True)
                    del mid_supra_parc

                    if "rh2refill" in locals():
                        #  Remove the voxels that are in the mid_parc and are in the rh2refill
                        ind = np.where((mid_parc.data != 0) & (rh2refill != 0))
                        rh2refill[ind] = 0

                    if "lh2refill" in locals():
                        #  Remove the voxels that are in the mid_parc and are in the lh2refill
                        ind = np.where((mid_parc.data != 0) & (lh2refill != 0))
                        lh2refill[ind] = 0

                    # mid_parc.save_parcellation(out_file= '/home/yaleman/mid_test.nii.gz', save_lut=True)

            # Detecting the number of regions
            if "lh_parc" in locals():
                nlh_subc = len(lh_parc.index)

            if "rh_parc" in locals():
                nrh_subc = len(rh_parc.index)

            if "mid_parc" in locals():
                nmid_subc = len(mid_parc.index)

            # if "WhiteMatter" in supra_names:
            # self.supra_dict[supra][supra][atlas_code]["mid"]["index"]

            date_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            if bool_ctx:

                # Atributes for the cortical parcellation
                atlas_names = self.parc_dict["Cortical"]["parcels"]

                proc_dict = self.parc_dict["Cortical"]["processing"]
                ctx_meth = proc_dict["method"]

                nctx_parc = len(
                    self.parc_dict["Cortical"]["processing"]["labels"]["lh"]
                )
                for c in np.arange(nctx_parc):

                    # Temporal header lines
                    glob_header_info_tmp = copy.deepcopy(glob_header_info)

                    ## -------- Cortical parcellation for the left hemisphere ---------------
                    # Creating the name for the output file
                    lh_in_parc = self.parc_dict["Cortical"]["processing"]["labels"][
                        "lh"
                    ][c]
                    at_name = [s for s in atlas_names if s in lh_in_parc]
                    lh_out_annot = os.path.join(
                        deriv_dir,
                        self.parc_dict["Cortical"]["deriv_surffold"],
                        path_cad,
                        "anat",
                        fullid + "_hemi-L" + "_" + "".join(at_name) + "_dseg.annot",
                    )

                    ## -------- Cortical parcellation for the right hemisphere ---------------
                    # Creating the name for the output file
                    rh_in_parc = self.parc_dict["Cortical"]["processing"]["labels"][
                        "rh"
                    ][c]
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
                            ref_id=self.parc_dict["Cortical"]["processing"][
                                "reference"
                            ],
                            hemi="lh",
                            fs_annot=lh_in_parc,
                            ind_annot=lh_out_annot,
                            cont_tech=cont_tech_freesurfer,
                            cont_image=cont_image_freesurfer,
                            force=force,
                        )

                        sub2proc.annot2ind(
                            ref_id=self.parc_dict["Cortical"]["processing"][
                                "reference"
                            ],
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

                    ## -------- Creating the volumetric parcellation ---------------
                    out_vol_dir = os.path.join(
                        deriv_dir,
                        self.parc_dict["Cortical"]["deriv_volfold"],
                        path_cad,
                        "anat",
                    )
                    if growwm is None:
                        growwm = ["0"]

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

                    for ngrow in np.arange(len(growwm)):
                        if growwm[ngrow] == "0":
                            out_vol_name = fullid + "_" + at_name + "_dseg.nii.gz"
                        else:

                            ent_dict = cltbids.str2entity(at_name)
                            if "desc" in ent_dict.keys():
                                if bool_mixwm:
                                    ent_dict["desc"] = (
                                        ent_dict["desc"]
                                        + "grow"
                                        + str(growwm[ngrow])
                                        + "mm+mixwm"
                                    )
                                else:
                                    ent_dict["desc"] = (
                                        ent_dict["desc"]
                                        + "grow"
                                        + str(growwm[ngrow])
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
                                        + str(growwm[ngrow])
                                        + "mm+mixwm_dseg.nii.gz"
                                    )
                                else:
                                    out_vol_name = (
                                        fullid
                                        + "_"
                                        + at_name
                                        + "_desc-grow"
                                        + str(growwm[ngrow])
                                        + "mm_dseg.nii.gz"
                                    )

                        sub2proc.surf2vol(
                            atlas=at_name,
                            out_vol=os.path.join(out_vol_dir, out_vol_name),
                            gm_grow=growwm[ngrow],
                            bool_mixwm=bool_mixwm,
                            force=False,
                            bool_native=True,
                            color_table=["tsv", "lut"],
                            cont_tech=cont_tech_freesurfer,
                            cont_image=cont_image_freesurfer,
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

                        # Creating the first part of the headers
                        part_header = [
                            "# $Id: {} {} \n".format(chim_parc_lut, date_time)
                        ]

                        part_header.append(
                            "# Corresponding parcellation: {} \n".format(chim_parc_file)
                        )

                        lut_header = part_header + glob_header_info_tmp
                        lut_header = lut_header + ["\n"]
                        lut_header.append(
                            "{:<4} {:<50} {:>3} {:>3} {:>3} {:>3}".format(
                                "#No.", "Label Name:", "R", "G", "B", "A"
                            )
                        )

                        if (
                            not os.path.isfile(chim_parc_file)
                            or not os.path.isfile(chim_parc_lut)
                            or not os.path.isfile(chim_parc_tsv)
                            or force
                        ):
                            # Creating the joined parcellation
                            ref_image = np.zeros_like(
                                t1_image.get_fdata(), dtype=np.int32
                            )
                            chim_parc = cltparc.Parcellation(
                                parc_file=ref_image, affine=affine
                            )

                            ctx_parc = cltparc.Parcellation(
                                parc_file=os.path.join(out_vol_dir, out_vol_name)
                            )
                            ctx_parc.remove_by_name(
                                names2remove=["unknown", "medialwall", "corpuscallosum"]
                            )

                            lh_ctx_parc = copy.deepcopy(ctx_parc)
                            rh_ctx_parc = copy.deepcopy(ctx_parc)

                            lh_ctx_parc.keep_by_name(names2look="ctx-lh-")
                            nlh_ctx = len(lh_ctx_parc.index)
                            rh_ctx_parc.keep_by_name(names2look="ctx-rh-")
                            nrh_ctx = len(rh_ctx_parc.index)

                            # Detect the global White Matter
                            brain_wm_parc = copy.deepcopy(ctx_parc)
                            brain_wm_parc.keep_by_code(
                                codes2keep=[2, 41, 5001, 5002, 7, 46]
                            )
                            ind = np.where(brain_wm_parc.data != 0)
                            brain_wm_parc.data[ind] = 1
                            brain_wm_parc.index = [1]
                            brain_wm_parc.name = ["wm-brain-whitematter"]
                            brain_wm_parc.color = ["#ffffff"]
                            brain_wm_parc.rearrange_parc(offset=2999)
                            brain_wm_parc.data[np.where(lh2refill)] = 3000
                            brain_wm_parc.data[np.where(rh2refill)] = 3000

                            # White Matter for the Right Hemisphere
                            tmp_rh = cltmisc.filter_by_substring(
                                ctx_parc.name, "wm-rh-"
                            )
                            if tmp_rh:
                                rh_wm_parc = copy.deepcopy(ctx_parc)
                                rh_wm_parc.keep_by_name(names2look=tmp_rh)
                                rh_wm_parc.rearrange_parc(offset=3000)

                            # White Matter for the Left Hemisphere
                            tmp_lh = cltmisc.filter_by_substring(
                                ctx_parc.name, "wm-lh-"
                            )
                            if tmp_lh:
                                lh_wm_parc = copy.deepcopy(ctx_parc)
                                lh_wm_parc.keep_by_name(names2look=tmp_lh)
                                lh_wm_parc.rearrange_parc(
                                    offset=3000 + nrh_ctx + nrh_subc
                                )

                            # Adding the right cortical parcellation to the final image
                            rh_ctx_parc.rearrange_parc()
                            chim_parc.add_parcellation(rh_ctx_parc, append=True)
                            del rh_ctx_parc

                            # Adding the right non-cortical parcellation to the final image
                            if "rh_parc" in locals():
                                chim_parc.add_parcellation(rh_parc, append=True)

                            # Adding the left cortical parcellation to the final image
                            lh_ctx_parc.rearrange_parc()
                            chim_parc.add_parcellation(lh_ctx_parc, append=True)
                            del lh_ctx_parc

                            # Adding the left non-cortical parcellation to the final image
                            if "lh_parc" in locals():
                                chim_parc.add_parcellation(lh_parc, append=True)

                            # Adding the regions that do not belong to any hemisphere to the final image
                            if "mid_parc" in locals():
                                chim_parc.add_parcellation(mid_parc, append=True)

                            # Adding the white matter to the final image
                            chim_parc.add_parcellation(brain_wm_parc, append=False)
                            del brain_wm_parc

                            if "rh_wm_parc" in locals():
                                chim_parc.add_parcellation(rh_wm_parc, append=False)
                                del rh_wm_parc

                            if "lh_wm_parc" in locals():
                                chim_parc.add_parcellation(lh_wm_parc, append=False)
                                del lh_wm_parc

                            # Adding the extra regions
                            if "extra_parc" in locals():
                                # Detecting if there is region overlap and removing it
                                tmp_extra = copy.deepcopy(extra_parc)
                                mask = np.logical_and(
                                    tmp_extra.data != 0, chim_parc.data != 0
                                )
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
                                save_lut=True,
                                save_tsv=True,
                            )
                            del chim_parc
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
                    part_header = ["# $Id: {} {} \n".format(chim_parc_lut, date_time)]
                    part_header.append(
                        "# Corresponding parcellation: {} \n".format(chim_parc_file)
                    )

                    lut_header = part_header + glob_header_info
                    lut_header = lut_header + ["\n"]
                    lut_header.append(
                        "{:<4} {:<50} {:>3} {:>3} {:>3} {:>3}".format(
                            "#No.", "Label Name:", "R", "G", "B", "A"
                        )
                    )

                    # Creating the joined parcellation
                    ref_image = np.zeros_like(t1_image.get_fdata(), dtype=np.int32)
                    chim_parc = cltparc.Parcellation(parc_file=ref_image, affine=affine)

                    # Adding the right non-cortical parcellation to the final image
                    if "rh_parc" in locals():
                        rh_parc.rearrange_parc()
                        chim_parc.add_parcellation(rh_parc, append=True)

                    # Adding the left non-cortical parcellation to the final image
                    if "lh_parc" in locals():
                        lh_parc.rearrange_parc()
                        chim_parc.add_parcellation(lh_parc, append=True)

                    # Adding the regions that do not belong to any hemisphere to the final image
                    if "mid_parc" in locals():
                        mid_parc.rearrange_parc()
                        chim_parc.add_parcellation(mid_parc, append=True)

                    # Saving the FINAL parcellation
                    chim_parc.save_parcellation(
                        out_file=chim_parc_file,
                        affine=affine,
                        headerlines=lut_header,
                        save_lut=True,
                        save_tsv=True,
                    )
                    del chim_parc


def _build_args_parser():

    class ColoredHelpFormatter(argparse.HelpFormatter):
        def __init__(self, prog):
            super().__init__(prog, max_help_position=52, width=100)
        
        def _format_action_invocation(self, action):
            # This method formats the option strings part (e.g. "--bidsdir PATH, -b PATH")
            if not action.option_strings:
                return super()._format_action_invocation(action)
            
            parts = []
            for option_string in action.option_strings:
                if option_string.startswith('--'):
                    colored_option = f"{bcolors.BOLD}{bcolors.OKBLUE}{option_string}{bcolors.ENDC}"
                elif option_string.startswith('-'):
                    colored_option = f"{bcolors.OKYELLOW}{option_string}{bcolors.ENDC}"
                else:
                    colored_option = option_string
                
                # Add metavar if present
                if action.metavar:
                    colored_option += f" {action.metavar}"
                parts.append(colored_option)
            
            return ', '.join(parts)
        
        def start_section(self, heading):
            if heading:
                if "Required" in heading:
                    heading = f"{bcolors.BOLD}{bcolors.OKGREEN}═══ {heading.upper()} ═══{bcolors.ENDC}"
                elif "Optional" in heading:
                    heading = f"{bcolors.BOLD}{bcolors.DARKCYAN}═══ {heading.upper()} ═══{bcolors.ENDC}"
                else:
                    heading = f"{bcolors.BOLD}{heading}{bcolors.ENDC}"
            super().start_section(heading)
    
    from argparse import ArgumentParser
    
    description = f"""
{bcolors.BOLD}{bcolors.HEADER}╔══════════════════════════════════════════════════════════════╗{bcolors.ENDC}
{bcolors.BOLD}{bcolors.HEADER}║                      CHIMERA TOOL                           ║{bcolors.ENDC}
{bcolors.BOLD}{bcolors.HEADER}║         Brain Parcellation Fusion Framework                 ║{bcolors.ENDC}
{bcolors.BOLD}{bcolors.HEADER}╚══════════════════════════════════════════════════════════════╝{bcolors.ENDC}

{bcolors.ITALIC}{bcolors.DARKWHITE}Generate custom brain parcellations by combining multiple atlases for different brain regions.{bcolors.ENDC}
"""

    epilog = f"""
{bcolors.BOLD}Examples:{bcolors.ENDC}
  chimera --regions                           # List available parcellations
  chimera -b /data -p HFMIIIFIFN -d /out -fr /fs  # Basic usage with 10-char code

{bcolors.BOLD}For more information:{bcolors.ENDC}
  Use --regions to see all available parcellation codes for each brain region.
"""

    p = argparse.ArgumentParser(
        formatter_class=ColoredHelpFormatter, 
        description=description,
        epilog=epilog
    )
    
    # Required arguments group (only the 3 truly required ones)
    requiredNamed = p.add_argument_group("Required Arguments")
    
    requiredNamed.add_argument(
        "--bidsdir",
        "-b",
        action="store",
        required=False,
        metavar="PATH",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}BIDS dataset directory{bcolors.ENDC}\n"
             f"Path to the Brain Imaging Data Structure (BIDS) dataset folder.\n",
    )
    
    requiredNamed.add_argument(
        "--parcodes",
        "-p",
        action="store",
        required=False,
        metavar="CODE",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}Parcellation code sequence{bcolors.ENDC}\n"
             f"10-character string identifying parcellation for each brain region:\n"
             f"  {bcolors.OKYELLOW}1) Cortex, 2) Basal ganglia, 3) Thalamus, 4) Amygdala, 5) Hippocampus{bcolors.ENDC}\n"
             f"  {bcolors.OKYELLOW}6) Hypothalamus, 7) Cerebellum, 8) Brainstem, 9) Gyral WM, 10) WM{bcolors.ENDC}\n\n"
             f"{bcolors.UNDERLINE}Example:{bcolors.ENDC} HFMIIIFIFN\n"
             f"Use {bcolors.OKGREEN}--regions{bcolors.ENDC} to see all available codes.\n",
    )
    
    requiredNamed.add_argument(
        "--subjids",
        "-ids",
        action="store",
        required=False,
        metavar="IDS",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}Subject identifiers{bcolors.ENDC}\n"
             f"Comma-separated subject IDs or path to text file with IDs.\n"
             f"{bcolors.ITALIC}Example file format:{bcolors.ENDC}\n"
             f"  sub-00001_ses-0001_run-2\n"
             f"  sub-00001_ses-0003_run-1\n",
        default=None,
    )
    
    # Optional arguments group  
    optionalNamed = p.add_argument_group("Optional Arguments")
    
    optionalNamed.add_argument(
        "--regions",
        "-r",
        action="store_true",
        required=False,
        help=f"{bcolors.BOLD}List available parcellations{bcolors.ENDC}\n"
             f"Display all parcellation options for each brain region.\n",
    )
    
    optionalNamed.add_argument(
        "--derivdir",
        "-d",
        action="store",
        required=False,
        metavar="PATH",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}Derivatives directory{bcolors.ENDC}\n"
             f"Output folder for results. Created if it doesn't exist.\n"
             f"If not specified, creates 'derivatives' inside BIDS directory.\n",
        default=None,
    )
    
    optionalNamed.add_argument(
        "--freesurferdir",
        "-fr",
        action="store",
        required=False,
        metavar="PATH",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}FreeSurfer subjects directory{bcolors.ENDC}\n"
             f"Path to FreeSurfer SUBJECTS_DIR. Created if it doesn't exist.\n",
        default=None,
    )
    
    optionalNamed.add_argument(
        "--scale",
        "-s",
        action="store",
        required=False,
        metavar="SCALE",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}Scale identifier{bcolors.ENDC}\n"
             f"Required for multi-resolution parcellations (e.g., Lausanne, Schaefer).\n"
             f"If not specified, generates all available scales.\n",
        default=None,
    )
    
    optionalNamed.add_argument(
        "--seg",
        "-e",
        action="store",
        required=False,
        metavar="SEG",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}Segmentation identifier{bcolors.ENDC}\n"
             f"Required when parcellations have multiple versions\n"
             f"(e.g., Schaefer: '7n' vs 'kong7n').\n",
        default=None,
    )
    
    optionalNamed.add_argument(
        "--nthreads",
        "-n",
        action="store",
        required=False,
        metavar="N",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}Number of parallel processes{bcolors.ENDC}\n"
             f"Number of subjects to process simultaneously (default: 1).\n",
        default=["1"],
    )
    
    optionalNamed.add_argument(
        "--growwm",
        "-g",
        action="store",
        required=False,
        metavar="MM",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}White matter growth{bcolors.ENDC}\n"
             f"Expand gray matter labels into white matter (in mm).\n"
             f"Multiple values can be comma-separated.\n",
        default=None,
    )
    
    optionalNamed.add_argument(
        "--config",
        "-c",
        action="store",
        required=False,
        metavar="FILE",
        type=str,
        nargs=1,
        help=f"{bcolors.BOLD}Pipeline configuration file{bcolors.ENDC}\n"
             f"Custom configuration file for advanced settings.\n",
        default=None,
    )
    
    optionalNamed.add_argument(
        "--mergectx",
        "-mctx",
        action="store_true",
        required=False,
        help=f"{bcolors.BOLD}Merge cortical regions{bcolors.ENDC}\n"
             f"Combine cortical gray matter and white matter regions.\n",
        default=False,
    )
    
    optionalNamed.add_argument(
        "--force",
        "-f",
        action="store_true",
        required=False,
        help=f"{bcolors.BOLD}Force overwrite{bcolors.ENDC}\n"
             f"Overwrite existing results without prompting.\n",
    )
    
    optionalNamed.add_argument(
        "--verbose",
        "-v",
        action="store",
        required=False,
        type=int,
        nargs=1,
        metavar="LEVEL",
        help=f"{bcolors.BOLD}Verbosity level{bcolors.ENDC}\n"
             f"0=quiet, 1=normal, 2=debug (default: 0).\n",
    )

    args = p.parse_args()

    global bids_dirs, supra_dict, deriv_dirs, fssubj_dirs, parcodes, pipe_json

    pipe_json = args.config

    if isinstance(args.config, list):
        pipe_json = args.config[0]

    if args.regions is True:
        if args.bidsdir is None and args.parcodes is None:
            print("\n")
            mess = "Available parcellations for each supra-region"
            print(
                "{}{}{}{}{}: ".format(
                    bcolors.BOLD, bcolors.PURPLE, mess, bcolors.ENDC, bcolors.ENDC
                )
            )
            _print_availab_parcels()
            sys.exit()

        elif args.bidsdir is None or args.parcodes is None:
            print("--bidsdir and --parcodes are REQUIRED arguments")
            sys.exit()

    bids_dirs = args.bidsdir[0].split(sep=",")
    # Remove possible empty elements
    bids_dirs = [x for x in bids_dirs if x]

    for bids_dir in bids_dirs:
        if not os.path.isdir(bids_dir):
            print("Please, supply a valid BIDs directory.")
            print(
                "The supplied BIDs directory does not exist: {}{}{}{}{}: is not supplied. ".format(
                    bcolors.BOLD, bcolors.OKRED, bids_dir, bcolors.ENDC, bcolors.ENDC
                )
            )
            p.print_help()
            sys.exit()

    if args.derivdir is None:
        print(
            "{}{}{}{}{}: is not supplied. ".format(
                bcolors.BOLD,
                bcolors.OKMAGENTA,
                "--derivdir",
                bcolors.ENDC,
                bcolors.ENDC,
            )
        )
        print(
            "The derivatives directory will be created in the corresponding BIDs directory."
        )
        deriv_dirs = []
        for bids_dir in bids_dirs:
            # print('derivatives_dir: {}'.format(os.path.join(bids_dir, 'derivatives')))
            print(
                "{}{}{}{}{}: {}{}{} ".format(
                    bcolors.BOLD,
                    bcolors.DARKCYAN,
                    "derivatives_dir",
                    bcolors.ENDC,
                    bcolors.ENDC,
                    bcolors.UNDERLINE,
                    os.path.join(bids_dir, "derivatives"),
                    bcolors.ENDC,
                )
            )

            deriv_dir = Path(os.path.join(bids_dir, "derivatives"))
            deriv_dir.mkdir(parents=True, exist_ok=True)

            deriv_dirs.append(str(deriv_dir))

    else:
        deriv_dirs = args.derivdir[0].split(sep=",")
        # Remove possible empty elements
        deriv_dirs = [x for x in deriv_dirs if x]

        if len(deriv_dirs) != len(bids_dirs):
            print(
                "The number of derivatives directories should be the same as the number of BIDs directories."
            )
            print(
                "The first derivatives directory will be the same for all BIDs directories:"
            )
            # print('derivatives_dir: {}'.format(deriv_dirs[0]))
            print(
                "{}{}{}{}{}: {}{}{} ".format(
                    bcolors.BOLD,
                    bcolors.DARKCYAN,
                    "derivatives_dir",
                    bcolors.ENDC,
                    bcolors.ENDC,
                    bcolors.UNDERLINE,
                    deriv_dirs[0],
                    bcolors.ENDC,
                )
            )

            deriv_dir = Path(deriv_dirs[0])

            # Create the folder if it does not exist
            deriv_dir.mkdir(parents=True, exist_ok=True)

            # Create a list of the same length as the number of BIDs directories
            deriv_dirs = [str(deriv_dir) for i in range(len(bids_dirs))]
        else:
            for deriv_dir in deriv_dirs:
                deriv_dir = Path(deriv_dir)

                # Create the folder if it does not exist
                deriv_dir.mkdir(parents=True, exist_ok=True)
    print("\n")
    if args.freesurferdir is None:

        print(
            "{}{}{}{}{}: is not supplied. ".format(
                bcolors.BOLD,
                bcolors.OKMAGENTA,
                "--freesurferdir",
                bcolors.ENDC,
                bcolors.ENDC,
            )
        )

        if "SUBJECTS_DIR" in os.environ:
            print(
                "The FreeSurfer subjects directory will be the same for all derivatives directories."
            )
            print("We will use the enviroment variable SUBJECTS_DIR.")
            print(
                "{}{}{}{}{}: {}{}{} ".format(
                    bcolors.BOLD,
                    bcolors.DARKCYAN,
                    "freesurfer_dir",
                    bcolors.ENDC,
                    bcolors.ENDC,
                    bcolors.UNDERLINE,
                    os.environ["SUBJECTS_DIR"],
                    bcolors.ENDC,
                )
            )

            fssubj_dir = Path(os.environ["SUBJECTS_DIR"])
            fssubj_dir.mkdir(parents=True, exist_ok=True)
            fssubj_dirs = [str(fssubj_dir) for i in range(len(deriv_dirs))]

        else:
            print(
                "The FreeSurfer subjects directory will be created in the following derivatives directory:"
            )
            fssubj_dirs = []
            for deriv_dir in deriv_dirs:
                print(
                    "freesurfer_dir: {}".format(os.path.join(deriv_dir, "freesurfer"))
                )
                print(
                    "{}{}{}{}{}: {}{}{} ".format(
                        bcolors.BOLD,
                        bcolors.DARKCYAN,
                        "freesurfer_dir",
                        bcolors.ENDC,
                        bcolors.ENDC,
                        bcolors.UNDERLINE,
                        os.path.join(deriv_dir, "freesurfer"),
                        bcolors.ENDC,
                    )
                )
                fssubj_dir = Path(os.path.join(deriv_dir, "freesurfer"))
                fssubj_dir.mkdir(parents=True, exist_ok=True)

                fssubj_dirs.append(str(fssubj_dir))
    else:
        fssubj_dirs = args.freesurferdir[0].split(sep=",")
        # Remove possible empty elements
        fssubj_dirs = [x for x in fssubj_dirs if x]

        if len(fssubj_dirs) != len(deriv_dirs):
            print(
                "The number of freesurfer directories should be the same as the number of derivatives directories."
            )
            print(
                "The FreeSurfer subjects directory  will be the same for all derivatives directories"
            )
            print(
                "{}{}{}{}{}: {}{}{} ".format(
                    bcolors.BOLD,
                    bcolors.DARKCYAN,
                    "freesurfer_dir",
                    bcolors.ENDC,
                    bcolors.ENDC,
                    bcolors.UNDERLINE,
                    fssubj_dirs[0],
                    bcolors.ENDC,
                )
            )

            fssubj_dir = Path(fssubj_dirs[0])

            # Create the folder if it does not exist
            fssubj_dir.mkdir(parents=True, exist_ok=True)

            # Create a list of the same length as the number of BIDs directories
            fssubj_dirs = [str(fssubj_dir) for i in range(len(deriv_dirs))]
        else:
            for fssubj_dir in fssubj_dirs:
                fssubj_dir = Path(fssubj_dir)

                # Create the folder if it does not exist
                fssubj_dir.mkdir(parents=True, exist_ok=True)

    parcodes = args.parcodes[0].split(sep=",")
    parc_dict, supra_dict = load_parcellations_info()
    supra_reg_names = list(parc_dict.keys())
    n_supra = len(supra_reg_names)

    # Remove empty elements
    parcodes = [x for x in parcodes if x]
    for i, parcode in enumerate(parcodes):
        if len(parcode) != n_supra:
            parcode = parcode.ljust(n_supra, "N")
            parcodes[i] = parcode

        # Checking if the code is correct
        bool_exit = False

        for ord, sp in enumerate(supra_reg_names):
            if parcode[ord] not in parc_dict[sp].keys() and parcode[ord] != "N":
                bool_exit = True
                print(
                    f"The parcellation code for the {sp} ({parcode[ord]}) is not correct."
                )
                _print_availab_parcels(sp)
                print(" ")
                print(
                    f"The {sp} structures will not me included in the final parcellation"
                )
                print(" ")

    return p


# simple progress indicator callback function
def progress_indicator(future):
    """
    A simple progress indicator for the concurrent futures
    :param future: future object

    """
    global lock, n_subj, n_comp, pb, pb1, chim_code
    # obtain the lock
    with lock:
        # update the counter
        n_comp += 1
        # report progress
        # print(f'{tasks_completed}/{n_subj} completed, {n_subj-tasks_completed} remain.')
        # pb.update(task_id=pb1, description= f'[red]Completed {n_comp}/{n_subj}', completed=n_subj)
        pb.update(
            task_id=pb1,
            description=f"[red]{chim_code}: Finished ({n_comp}/{n_subj})",
            completed=n_comp,
        )


def chimera_parcellation(
    bids_dir: str,
    deriv_dir: str,
    fssubj_dir: str,
    code_dict: dict,
    t1s2run_file: str = None,
    growwm: list = ["0"],
    mixwm: bool = False,
    nthreads: int = 1,
):
    """
    Preparing chimera to build the parcellations.

    Parameters
    ----------
    bids_dir : str
        The directory with the input dataset formatted according to the BIDS standard.

    deriv_dir : str
        The directory where the output files are stored.

    code_dict : dict
        Dictionary containing 3 main keys: code, scale and seg
        -- "code" can be either a string or a list of strings.
        -- "scale" can be either a string or a list of strings.
        -- "seg" can be either a string or a list of strings.

    t1s2run_file : str
        File containing the list of T1w images to be processed.

    growwm : list
        List of values, in mm, to grow the GM regions inside the WM.

    nthreads : int
        Number of threads to be used.

    """

    # Declaring global variables
    global pipe_json, pipe_dict, pb, pb1, n_subj, n_comp, lock, chim_code

    ######## -- Reading the configuration dictionary  ------------ #
    pipe_dict = _pipeline_info(pipe_json=pipe_json)

    if t1s2run_file is None:
        # Selecting all the T1w images for each BIDS directory
        layout = BIDSLayout(bids_dir, validate=False, derivatives=False)
        t1s = layout.get(
            extension=["nii.gz", "nii"], suffix="T1w", return_type="filename"
        )
    else:
        if os.path.exists(t1s2run_file):
            with open(t1s2run_file) as file:
                t1s2run = [line.rstrip() for line in file]
        else:

            # If the string contains a file separator, we assume that the string is a list of files
            if os.path.sep not in t1s2run_file:
                t1s2run = t1s2run_file.split(",")
            else:
                raise ValueError(
                    "Please, provide a valid file containing the list of T1w images to be processed."
                )

        t1s = []
        for id in t1s2run:

            if not os.path.isfile(id):
                id_ent = cltbids.str2entity(id)
                if "ses" in id_ent.keys():
                    path_cad = os.path.join(
                        bids_dir, "sub-" + id_ent["sub"], "ses-" + id_ent["ses"], "anat"
                    )
                else:
                    path_cad = os.path.join(bids_dir, "sub-" + id_ent["sub"], "anat")

                if "suffix" not in id_ent.keys():
                    id_ent["suffix"] = "T1w"

                if "extension" not in id_ent.keys():
                    id_ent["extension"] = "nii.gz"

                t1_temp = os.path.join(path_cad, cltbids.entity2str(id_ent))
                if os.path.isfile(t1_temp):
                    t1s.append(t1_temp)

    chim_codes = code_dict["code"]

    n_parc = len(chim_codes)
    n_subj = len(t1s)

    with Progress() as pb:
        pb2 = pb.add_task("[green]Parcellation: ", total=n_parc)

        # Loop around each parcellation
        for p, chim_code in enumerate(chim_codes):

            # Creating and configuring the Chimera object
            chim_obj = Chimera(
                parc_code=chim_code, scale=code_dict["scale"], seg=code_dict["seg"]
            )

            # Creating the color table
            chim_obj.create_table()

            # Configuring and downloading the templates
            chim_obj.prepare_templates(fssubj_dir=fssubj_dir)

            # create a lock for the counter
            lock = Lock()

            # Completed subjects
            n_comp = 0
            failed = []

            pb.update(
                task_id=pb2,
                description=f"[green]Parcellation: {chim_code} ({p+1}/{n_parc})",
                completed=p + 1,
            )

            # print("Parcellation: % d"% (p+1), "of % d"% (n_parc))
            if nthreads == 1:

                pb1 = pb.add_task(
                    f"[red]Processing: Subject ({1}/{n_subj}) ", total=n_subj
                )
                for i, t1 in enumerate(t1s):
                    # ent_dict = layout.parse_file_entities(t1)

                    t1_name = os.path.basename(t1)
                    temp = t1_name.split("_")
                    full_id = "_".join(temp[:-1])
                    pb.update(
                        task_id=pb1,
                        description=f"[red]{chim_code}: {full_id} ({i+1}/{n_subj})",
                        completed=i + 1,
                    )

                    chim_obj.build_parcellation(
                        t1, bids_dir, deriv_dir, fssubj_dir, growwm, mixwm
                    )

            else:
                start_time = time.perf_counter()

                # create a progress bar for the subjects
                pb1 = pb.add_task(
                    f"[red]Processing: Subject ({1}/{n_subj}) ", total=n_subj
                )

                # Adjusting the number of threads to the number of subjects
                if n_subj < nthreads:
                    nthreads = n_subj

                # start the thread pool
                with ThreadPoolExecutor(nthreads) as executor:
                    # send in the tasks
                    # futures = [executor.submit(_build_parcellation, t1s[i],
                    # bids_dir, deriv_dir, parccode, growwm, mixwm) for i in range(n_subj)]

                    futures = [
                        executor.submit(
                            chim_obj.build_parcellation,
                            t1s[i],
                            bids_dir,
                            deriv_dir,
                            fssubj_dir,
                            growwm,
                            mixwm,
                        )
                        for i in range(n_subj)
                    ]

                    # register the progress indicator callback
                    for future in futures:
                        future.add_done_callback(progress_indicator)
                    # wait for all tasks to complete


def main():
    # 0. Handle inputs
    parser = _build_args_parser()
    args = parser.parse_args()

    print(args)
    if args.verbose is not None:
        v = np.int(args.verbose[0])
    else:
        v = 0
        print("- Verbose set to 0\n")
    if v:
        print("\nInputs\n")
    #

    global bids_dirs, deriv_dirs, fssubj_dirs, parcodes

    if args.scale is not None:
        scale_id = args.scale[0].split(sep=",")

        # Remove empty elements
        scale_id = [x for x in scale_id if x]

    else:
        scale_id = None

    if args.seg is not None:
        seg_id = args.seg[0].split(sep=",")

        # Remove empty elements
        seg_id = [x for x in seg_id if x]
    else:
        seg_id = None

    # Create dictionary with the code, scale and segmentation
    code_dict = {"code": parcodes, "scale": scale_id, "seg": seg_id}

    if args.subjids is not None:
        t1s2run_file = args.subjids[0]
    else:
        t1s2run_file = None

    if args.growwm is not None:
        growwm = args.growwm[0].split(sep=",")

        # Remove empty elements
        growwm = [x for x in growwm if x]
    else:
        growwm = ["0"]

    mixwm = args.mergectx

    # Detecting the number of cores to be used
    ncores = os.cpu_count()
    nthreads = int(args.nthreads[0])

    if nthreads > ncores:
        if ncores > 3:
            nthreads = ncores - 2
        else:
            nthreads = 2

    for i, bids_dir in enumerate(bids_dirs):

        deriv_dir = deriv_dirs[i]
        fssubj_dir = fssubj_dirs[i]
        chimera_parcellation(
            bids_dir,
            deriv_dir,
            fssubj_dir,
            code_dict,
            t1s2run_file,
            growwm,
            mixwm,
            nthreads,
        )


if __name__ == "__main__":
    main()
