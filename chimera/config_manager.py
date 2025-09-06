import os
import json
import numpy as np
import pandas as pd
from glob import glob

from pathlib import Path
import templateflow.api as tflow


# Loading the JSON file containing the available parcellations
def load_parcellations_info(parc_json: str = None, supra_folder: str = None):
    """
    Load the JSON file containing the available parcellations

    Parameters:
    ----------
    parc_json : str
        JSON file containing the parcellation dictionary.

    supra_folder : str
        Folder containing the supraregions TSV files.

    Returns:
    --------
    parc_dict : dict
        Dictionary containing the parcellation information

    supra_dict : dict
        Dictionary containing the supraregions information

    """
    chim_dir = os.path.dirname(os.path.abspath(__file__))

    # Get the absolute of this file
    if parc_json is None:
        parc_json = os.path.join(chim_dir, "config", "supraregions_dictionary.json")
    else:
        if not os.path.isfile(parc_json):
            raise ValueError(
                "Please, provide a valid JSON file containing the parcellation dictionary."
            )

    with open(parc_json) as f:
        parc_dict = json.load(f)

    # Read all the tsv files
    if supra_folder is None:
        spr_files = glob(os.path.join(chim_dir, "config", "supraregions", "*.tsv"))
    else:
        if os.path.isdir(supra_folder):
            spr_files = glob(os.path.join(supra_folder, "*.tsv"))
        else:
            raise ValueError(
                "Please, provide a valid folder containing the supraregions TSV files."
            )

    supra_dict = {}
    for spr in spr_files:
        spr_name = os.path.basename(spr).split(".")[0]
        temp_df = pd.read_csv(spr, sep="\t")
        sp_ids = temp_df["supraregion"].unique().tolist()

        # Create a dictionary for each supraregion
        supra_dict[spr_name] = {}

        sint_dict = {}  # Create a dictionary for each supraregion
        for sid in sp_ids:

            # Create a sub dataframe for each supraregion
            sub_df = temp_df.loc[temp_df["supraregion"] == sid]

            # For each supraregion id get the methods used to parcellate it
            st_methods = temp_df.loc[temp_df["supraregion"] == sid, "method"].tolist()
            st_methods = np.unique(st_methods).tolist()

            method_dict = {}  # Create a dictionary for each method
            for mid in st_methods:

                # Create a sub dataframe for each method
                sub_df2 = sub_df.loc[sub_df["method"] == mid]

                # Get the hemispheres
                st_hemi = sub_df2["hemi"].tolist()

                # Get the unique hemispheres
                st_hemi = np.unique(st_hemi).tolist()

                hemi_dict = {}  # Create a dictionary for each hemisphere
                for hemi in st_hemi:

                    # Create a sub dataframe for each hemisphere
                    sub_df3 = sub_df2.loc[sub_df2["hemi"] == hemi]

                    # Get the indexes of the regions
                    indexes = sub_df3["index"].tolist()

                    # Get the names of the regions
                    names = sub_df3["name"].tolist()

                    # Get the colors of the regions
                    colors = sub_df3["color"].tolist()

                    # Create a dictionary for each hemisphere
                    temp_dict = {"index": indexes, "name": names, "color": colors}

                    # Add the dictionary to the hemi_dict
                    hemi_dict[hemi] = temp_dict

                # Add the dictionary to the method_dict
                method_dict[mid] = hemi_dict

            # Add the dictionary to the supra_dict
            sint_dict[sid] = method_dict

        # Add the dictionary to the supra_dict
        supra_dict[spr_name] = sint_dict

    return parc_dict, supra_dict


# Loading the JSON file containing the available parcellations
def _pipeline_info(pipe_json: str = None):
    """
    Load the JSON file containing the pipeline configuration.

    Parameters:
    ----------
    pipe_json : str
        JSON file containing the pipeline configuration dictionary.

    Returns:
    --------
    pipe_dict : dict
        Dictionary containing the pipeline information

    """
    cwd = os.path.dirname(os.path.abspath(__file__))

    # Get the absolute of this file
    if pipe_json is None:

        pipe_json = os.path.join(cwd, "config", "pipe_config.json")
    else:
        if not os.path.isfile(pipe_json):
            raise ValueError(
                "Please, provide a valid JSON file containing the pipeline configuration dictionary."
            )

    with open(pipe_json, encoding="utf-8") as f:
        pipe_dict = json.load(f)

    return pipe_dict


def _set_templateflow_home(tflow_home: str = "local"):
    """
    Setting up the templateflow home directory.

    Parameters:
    ----------
    tflow_home : str or Path
        Templatefow home directory.

    Returns:
    --------
    updated_tflow_home : str
        Updated templateflow home directory.

    """

    orig_tflow_home = os.getenv("TEMPLATEFLOW_HOME")
    if tflow_home != "local":
        tflow_home = Path(tflow_home)
        tflow_home.mkdir(parents=True, exist_ok=True)

        if orig_tflow_home is None:
            if tflow_home.is_dir():
                updated_tflow_home = str(tflow_home)
                os.environ["TEMPLATEFLOW_HOME"] = updated_tflow_home
                tflow.update()
        else:
            if not os.path.samefile(orig_tflow_home, tflow_home):
                if tflow_home.is_dir():
                    updated_tflow_home = str(tflow_home)
                    os.environ["TEMPLATEFLOW_HOME"] = updated_tflow_home
                    tflow.update()
                else:
                    updated_tflow_home = orig_tflow_home
    else:
        updated_tflow_home = orig_tflow_home

    return updated_tflow_home
