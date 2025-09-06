import copy
import os
from glob import glob
import numpy as np
import pandas as pd

from .config_manager import load_parcellations_info

from clabtoolkit.misctools import bcolors
import clabtoolkit.parcellationtools as cltparc


def _mix_side_prop(st_dict: dict, boolsort: bool = True):
    """
    Method to mix all the properties for a specific supra-region and a specific method.

    Parameters:
    ----------
    st_dict : dict
        Dictionary containing the properties separated by sides.

    Returns:
    --------
    pipe_dict : dict
        Dictionary containing the pipeline information

    """

    all_index = []
    all_name = []
    all_color = []
    for side in st_dict.keys():
        all_index = all_index + st_dict[side]["index"]
        all_name = all_name + st_dict[side]["name"]
        all_color = all_color + st_dict[side]["color"]

    if boolsort:
        # Sort the all_index and apply the order to all_name and all_color
        sort_index = np.argsort(all_index)
        all_index = [all_index[i] for i in sort_index]
        all_name = [all_name[i] for i in sort_index]
        all_color = [all_color[i] for i in sort_index]

    return all_index, all_name, all_color


def create_extra_regions_parc(aparc: str, offset: int = 5000):
    """
    Create a parcellation object containing the extra regions. These parcellations
    are included in the Aparc+aseg image and include the regions stored in the Axiliary*.tsv files.
    These files are stored in the config/supraregions folder.

    Parameters:
    ----------
    aparc : str
        Aparc+aseg image obrained by FreeSurfer.

    Returns:
    --------
    extra_parc : Parcellation object
        Parcellation object containing the extra regions.

    """
    chim_dir = os.path.dirname(os.path.abspath(__file__))

    # Reading the auxiliary tsv file
    extra_tsv = glob(os.path.join(chim_dir, "config", "supraregions", "Auxiliary*.tsv"))

    # Reading the tsv file
    if extra_tsv:

        # Do a loop and read and concatenate all the files
        for ind, et in enumerate(extra_tsv):
            temp_df = pd.read_csv(et, sep="\t")
            # Append pandas dataframes
            if ind == 0:
                extra_df = temp_df
            else:
                extra_df = extra_df.append(temp_df, ignore_index=True)

    else:
        raise ValueError("Please, provide a valid auxiliary TSV file.")

    # Reading the Aparc+aseg image
    aparc_parc = cltparc.Parcellation(parc_file=aparc)

    # Create a dataframe for the left hemisphere
    lh_df = extra_df.loc[extra_df["hemi"] == "lh"]

    # Sort the dataframe according to the name
    lh_df = lh_df.sort_values(by="name")

    # Create a dataframe for the right hemisphere
    rh_df = extra_df.loc[extra_df["hemi"] == "rh"]

    # Sort the dataframe according to the name
    rh_df = rh_df.sort_values(by="name")

    # Create a dictionary for the structures without hemispheres
    mid_df = extra_df.loc[extra_df["hemi"] == "mid"]

    # Sort the dataframe according to the name
    mid_df = mid_df.sort_values(by="name")

    # Create the parcellation for the left hemisphere
    lh_tmp_parc = copy.deepcopy(aparc_parc)
    lh_tmp_parc.keep_by_code(codes2keep=lh_df["index"].tolist())
    lh_tmp_parc.index = lh_df["index"].tolist()
    lh_tmp_parc.name = lh_df["name"].tolist()
    lh_tmp_parc.color = lh_df["color"].tolist()
    lh_tmp_parc.adjust_values()
    lh_tmp_parc.rearrange_parc()

    # Create the parcellation for the right hemisphere
    rh_tmp_parc = copy.deepcopy(aparc_parc)
    rh_tmp_parc.keep_by_code(codes2keep=rh_df["index"].tolist())
    rh_tmp_parc.index = rh_df["index"].tolist()
    rh_tmp_parc.name = rh_df["name"].tolist()
    rh_tmp_parc.color = rh_df["color"].tolist()
    rh_tmp_parc.adjust_values()
    rh_tmp_parc.rearrange_parc()

    # Create the parcellation for the structures without hemispheres
    mid_tmp_parc = copy.deepcopy(aparc_parc)
    mid_tmp_parc.keep_by_code(codes2keep=mid_df["index"].tolist())
    mid_tmp_parc.index = mid_df["index"].tolist()
    mid_tmp_parc.name = mid_df["name"].tolist()
    mid_tmp_parc.color = mid_df["color"].tolist()
    mid_tmp_parc.adjust_values()
    mid_tmp_parc.rearrange_parc()

    # Unify the parcellations
    rh_tmp_parc.add_parcellation(lh_tmp_parc, append=True)
    rh_tmp_parc.add_parcellation(mid_tmp_parc, append=True)
    rh_tmp_parc.rearrange_parc(offset=offset)

    return rh_tmp_parc


def _print_availab_parcels(reg_name=None):
    """
    Print the available parcellations for each supra-region.

    Parameters:
    ----------
    reg_name : str
        Supra-region name. Default is None.

    Returns:
    --------

    """

    data, supra_dict = load_parcellations_info()

    if reg_name is None:
        supra_keys = data.keys()
        parc_help = ""
        for sup in supra_keys:
            parc_opts = data[sup]
            parc_help = '{} "{}:\n"'.format(parc_help, sup)
            print(
                "{}{}{}{}{}: ".format(
                    bcolors.BOLD, bcolors.DARKCYAN, sup, bcolors.ENDC, bcolors.ENDC
                )
            )

            for opts in parc_opts:
                desc = data[sup][opts]["name"]
                cita = data[sup][opts]["citation"]
                parc_help = '{} "{}: {} {}\n"'.format(parc_help, opts, desc, cita)
                print(
                    "{}     {}{}: {}{}{}{}{} {}{}{}".format(
                        bcolors.OKGREEN,
                        opts,
                        bcolors.ENDC,
                        bcolors.ITALIC,
                        bcolors.DARKWHITE,
                        desc,
                        bcolors.ENDC,
                        bcolors.ENDC,
                        bcolors.OKYELLOW,
                        cita,
                        bcolors.ENDC,
                    )
                )
            print("")
    else:
        parc_opts = data[reg_name]
        print("\n")
        print(
            "{}{}{}{}{}: ".format(
                bcolors.BOLD, bcolors.DARKCYAN, reg_name, bcolors.ENDC, bcolors.ENDC
            )
        )

        for opts in parc_opts:
            desc = data[reg_name][opts]["name"]
            cita = data[reg_name][opts]["citation"]
            print(
                "{}     {}{}: {}{}{}{}{} {}{}{}".format(
                    bcolors.OKGREEN,
                    opts,
                    bcolors.ENDC,
                    bcolors.ITALIC,
                    bcolors.DARKWHITE,
                    desc,
                    bcolors.ENDC,
                    bcolors.ENDC,
                    bcolors.OKYELLOW,
                    cita,
                    bcolors.ENDC,
                )
            )
        print("")
