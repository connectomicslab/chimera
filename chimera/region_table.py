#!/usr/bin/env python3
"""Region-table construction and export for the Chimera parcellation.

This module decomposes the original ``Chimera.create_table``,
``Chimera.export_table`` and ``Chimera.build_lut_header`` methods into named
functions.  Behaviour is preserved exactly; ``pipe_dict`` is threaded explicitly
instead of being read from a module-level global.
"""

import json
import os
from glob import glob
from os import PathLike
from pathlib import Path

import clabtoolkit.freesurfertools as cltfree
import clabtoolkit.misctools as cltmisc
import clabtoolkit.parcellationtools as cltparc
import numpy as np
import pandas as pd
from templateflow import api as tflow

from .config_manager import _pipeline_info, _set_templateflow_home


def build_region_table(
    chim,
    wm_index_offset: int = 3000,
    reg2rem: list | str | None = None,
    pipe_dict: dict = None,
):
    """Create the table of regions produced by the Chimera parcellation.

    The resulting tables are stored on ``chim.regtable`` exactly as the original
    ``Chimera.create_table`` method did.

    Parameters
    ----------
    chim : Chimera
        The Chimera object.

    wm_index_offset : int
        Index offset for the white matter parcellation (default = 3000).

    reg2rem : list
        List of regions to remove from the parcellation.

    pipe_dict : dict
        Pipeline configuration dictionary.  If ``None`` it is resolved with
        :func:`_pipeline_info`.
    """

    # Use a fresh list rather than a mutable default argument
    if reg2rem is None:
        reg2rem = ["unknown", "medialwall", "corpuscallosum"]

    if pipe_dict is None:
        pipe_dict = _pipeline_info()

    # Detecting the base directory
    cwd = os.path.dirname(os.path.abspath(__file__))
    chim_dir = os.path.dirname(cwd)

    # Reading the names of the supra-regions
    supra_names = list(chim.parc_dict.keys())

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
        atlas_code = chim.parc_dict[supra]["code"]
        atlas_str = chim.parc_dict[supra]["atlas"]
        atlas_desc = chim.parc_dict[supra]["description"]
        atlas_cita = chim.parc_dict[supra]["citation"]
        atlas_src = chim.parc_dict[supra]["source"]
        atlas_ref = chim.parc_dict[supra]["reference"]

        if supra == "Cortical":
            ctx_parc_lh, ctx_parc_rh = _load_cortical_table_files(
                chim, supra, atlas_src, atlas_str, atlas_ref, chim_dir, pipe_dict
            )

        else:

            desc_noctx.append(atlas_desc)
            atlas_type = chim.parc_dict[supra]["type"]
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

            if supra in chim.supra_dict.keys():
                meth_dict = chim.parc_dict[supra]
                st_dict = chim.supra_dict[supra][supra][meth_dict["code"]]
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

    # Removing the regions listed in reg2rem from the non-cortical structures
    _remove_regions_inplace(rh_noctx_names, rh_noctx_codes, rh_noctx_colors, reg2rem)
    _remove_regions_inplace(lh_noctx_names, lh_noctx_codes, lh_noctx_colors, reg2rem)
    _remove_regions_inplace(bs_noctx_names, bs_noctx_codes, bs_noctx_colors, reg2rem)

    # Creating the list of dataframes for the different parcellations
    if len(ctx_parc_lh) == 0:
        parc_id_list, desc_list, tab_list = _build_tables_without_cortex(
            chim,
            rh_noctx_names,
            lh_noctx_names,
            bs_noctx_names,
            rh_noctx_colors,
            lh_noctx_colors,
            bs_noctx_colors,
            desc_noctx,
        )
    else:
        parc_id_list, desc_list, tab_list = _build_tables_with_cortex(
            chim,
            ctx_parc_lh,
            ctx_parc_rh,
            reg2rem,
            wm_index_offset,
            desc_noctx,
            lh_noctx_names,
            rh_noctx_names,
            bs_noctx_names,
            lh_noctx_colors,
            rh_noctx_colors,
            bs_noctx_colors,
        )

    # Add the tab_list as an attribute of the class
    chim.regtable = {"parc_id": parc_id_list, "desc": desc_list, "table": tab_list}


def _remove_regions_inplace(names, codes, colors, reg2rem):
    """Remove, in place, the regions whose names match ``reg2rem``.

    Pop in descending order so earlier deletions do not shift the remaining
    indexes.
    """
    if names:
        indexes = cltmisc.get_indexes_by_substring(names, reg2rem).tolist()
        for i in sorted(indexes, reverse=True):
            names.pop(i)
            codes.pop(i)
            colors.pop(i)


def _load_cortical_table_files(
    chim, supra, atlas_src, atlas_str, atlas_ref, chim_dir, pipe_dict
):
    """Locate the cortical parcellation files used to build the region table."""

    # Atributes for the cortical parcellation
    atlas_type = chim.parc_dict[supra]["type"]
    atlas_names = chim.parc_dict[supra]["parcels"]

    ctx_parc_lh = []
    ctx_parc_rh = []

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
        else:
            raise ValueError(
                "The type of the cortical atlas is not valid. Please supply a valid type (annot or gcs)."
            )

        ctx_parc_lh = glob(os.path.join(atlas_dir, "*-L_*" + atlas_ext))
        ctx_parc_rh = glob(os.path.join(atlas_dir, "*-R_*" + atlas_ext))

        # Filtering for selecting the correct cortical parcellation
        ctx_parc_lh = cltmisc.filter_by_substring(
            ctx_parc_lh, or_filter=atlas_names, bool_case=False
        )
        ctx_parc_rh = cltmisc.filter_by_substring(
            ctx_parc_rh, or_filter=atlas_names, bool_case=False
        )

    return ctx_parc_lh, ctx_parc_rh


def _build_tables_without_cortex(
    chim,
    rh_noctx_names,
    lh_noctx_names,
    bs_noctx_names,
    rh_noctx_colors,
    lh_noctx_colors,
    bs_noctx_colors,
    desc_noctx,
):
    """Build the region table when no cortical parcellation is available."""

    tab_list = []
    desc_list = []
    parc_id = "atlas-chimera" + chim.parc_code
    parc_id_list = []

    all_names = rh_noctx_names + lh_noctx_names + bs_noctx_names
    all_colors = rh_noctx_colors + lh_noctx_colors + bs_noctx_colors
    index = np.arange(1, len(all_names) + 1).tolist()
    tab_df = pd.DataFrame({"index": index, "name": all_names, "color": all_colors})
    tab_list.append(tab_df)

    gen_desc = ["# Parcellation code: " + chim.parc_code]
    gen_desc.append(desc_noctx)

    parc_id_list.append(parc_id)

    return parc_id_list, desc_list, tab_list


def _build_tables_with_cortex(
    chim,
    ctx_parc_lh,
    ctx_parc_rh,
    reg2rem,
    wm_index_offset,
    desc_noctx,
    lh_noctx_names,
    rh_noctx_names,
    bs_noctx_names,
    lh_noctx_colors,
    rh_noctx_colors,
    bs_noctx_colors,
):
    """Build the region table when a cortical parcellation is available."""

    tab_list = []
    desc_list = []
    parc_id = "atlas-chimera" + chim.parc_code
    parc_id_list = []

    for i, parc_file in enumerate(ctx_parc_lh):

        gen_desc = ["# Parcellation code: " + chim.parc_code]

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

        gen_desc.append(chim.parc_dict["Cortical"]["description"])
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

        ## Removing elements from the table according to their name for both.
        ## Pop in descending order so earlier deletions do not shift indexes.
        indexes = cltmisc.get_indexes_by_substring(lh_ctx_name, reg2rem).tolist()
        for i in sorted(indexes, reverse=True):
            lh_ctx_name.pop(i)
            lh_ctx_color.pop(i)

        indexes = cltmisc.get_indexes_by_substring(rh_ctx_name, reg2rem).tolist()
        for i in sorted(indexes, reverse=True):
            rh_ctx_name.pop(i)
            rh_ctx_color.pop(i)

        # Concatenating the lists
        if "GyralWM" in chim.parc_dict.keys():
            gen_desc.append(chim.parc_dict["GyralWM"]["description"])

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
            lh_ctx_name + lh_noctx_names + bs_noctx_names + wm_rh_name + wm_lh_name
        )
        lh_all_indexes = np.arange(
            1, len(lh_ctx_name + lh_noctx_names + bs_noctx_names) + 1
        ) + np.max(rh_all_indexes)
        lh_all_indexes = lh_all_indexes.tolist()

        rh_all_colors = rh_ctx_color + rh_noctx_colors
        lh_all_colors = (
            lh_ctx_color + lh_noctx_colors + bs_noctx_colors + wm_rh_color + wm_lh_color
        )

        # Concatenating the hemispheres
        all_names = rh_all_names + lh_all_names
        all_colors = rh_all_colors + lh_all_colors
        all_indexes = rh_all_indexes + lh_all_indexes + wm_rh_indexes + wm_lh_indexes

        # Generating a dataframe
        tab_df = pd.DataFrame(
            {"index": all_indexes, "name": all_names, "color": all_colors}
        )
        tab_list.append(tab_df)
        desc_list.append(gen_desc)
        parc_id_list.append(parc_id)

    return parc_id_list, desc_list, tab_list


def export_table(chim, out_basename: str = None, format: list | str = "tsv"):
    """Export the table of regions to a TSV or a LUT file.

    Parameters
    ----------
    chim : Chimera
        The Chimera object holding the ``regtable`` attribute.

    out_basename : str
        Output basename for the TSV or LUT file.

    format : str or list
        Format of the output file. It can be 'tsv', 'lut' or ['tsv, lut'].
    """

    if out_basename is None:
        # Exit if the output basename is not provided
        raise ValueError("Please provide an output basename for the TSV or LUT file.")

    out_name = os.path.basename(out_basename)

    out_dir = os.path.dirname(out_basename)
    out_dir = Path(out_dir)

    # Create the output directory if it does not exist
    out_dir.mkdir(parents=True, exist_ok=True)

    parc_ids = chim.regtable["parc_id"]
    parc_desc = chim.regtable["desc"]
    parc_tables = chim.regtable["table"]

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


def build_lut_header(chim):
    """Build the header of the LUT file.

    Parameters
    ----------
    chim : Chimera
        The Chimera object.

    Returns
    -------
    list of str
        The header lines for the LUT file.
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
    chim_code = chim.parc_code

    # Creating the header lines
    headerlines = [f" # Chimera parcellation code: {chim.parc_code}"]
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
                f"    # {i + 1}. The parcellation code {chim_code[i]} is not present in the dictionary for the supra-region {supra}."
            )

    return headerlines
