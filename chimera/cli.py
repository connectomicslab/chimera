#!/usr/bin/env python3
"""Command-line interface and multi-subject orchestration for Chimera.

This module gathers the argument parser, the progress-indicator callback, the
``chimera_parcellation`` driver and ``main`` that previously lived at the bottom
of :mod:`chimera.chimera`.  It drives the refactored :class:`chimera.core.Chimera`.
"""

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Lock

import clabtoolkit.bidstools as cltbids
from bids import BIDSLayout
from clabtoolkit.colorstools import bcolors
from rich.progress import Progress, TaskID

from .config_manager import _pipeline_info, load_parcellations_info
from .core import Chimera
from .parcellation import _print_availab_parcels

# State shared between _build_args_parser, chimera_parcellation,
# progress_indicator and main via `global`. These are annotations only: they
# bind no value, so the runtime behaviour is unchanged and the names still come
# into existence when the functions that declare them global assign them. They
# exist so a reader (and a type checker) can see the shared state in one place
# instead of inferring it from scattered `global` statements.
bids_dirs: list[str]
deriv_dirs: list[str]
fssubj_dirs: list[str]
parcodes: list[str]
supra_dict: dict
pipe_json: str
pipe_dict: dict
lock: Lock
n_subj: int
n_comp: int
pb: Progress
pb1: TaskID
chim_code: str


def _build_args_parser():

    class ColoredHelpFormatter(argparse.RawDescriptionHelpFormatter):
        def __init__(self, prog):
            super().__init__(prog, max_help_position=52, width=100)

        def _format_action_invocation(self, action):
            # This method formats the option strings part (e.g. "--bidsdir PATH, -b PATH")
            if not action.option_strings:
                return super()._format_action_invocation(action)

            parts = []
            for option_string in action.option_strings:
                if option_string.startswith("--"):
                    colored_option = (
                        f"{bcolors.BOLD}{bcolors.OKBLUE}{option_string}{bcolors.ENDC}"
                    )
                elif option_string.startswith("-"):
                    colored_option = f"{bcolors.OKYELLOW}{option_string}{bcolors.ENDC}"
                else:
                    colored_option = option_string

                # Add metavar if present
                if action.metavar:
                    colored_option += f" {action.metavar}"
                parts.append(colored_option)

            return ", ".join(parts)

        def start_section(self, heading):
            if heading:
                if "Required" in heading:
                    heading = f"{bcolors.BOLD}{bcolors.OKGREEN}═══ {heading.upper()} ═══{bcolors.ENDC}"
                elif "Optional" in heading:
                    heading = f"{bcolors.BOLD}{bcolors.DARKCYAN}═══ {heading.upper()} ═══{bcolors.ENDC}"
                else:
                    heading = f"{bcolors.BOLD}{heading}{bcolors.ENDC}"
            super().start_section(heading)


    description = f"""
{bcolors.BOLD}{bcolors.HEADER}╔══════════════════════════════════════════════════════════════╗{bcolors.ENDC}
{bcolors.BOLD}{bcolors.HEADER}║                      CHIMERA TOOL                           ║{bcolors.ENDC}
{bcolors.BOLD}{bcolors.HEADER}║         Brain Parcellation Fusion Framework                 ║{bcolors.ENDC}
{bcolors.BOLD}{bcolors.HEADER}╚══════════════════════════════════════════════════════════════╝{bcolors.ENDC}

{bcolors.ITALIC}{bcolors.DARKWHITE}Generate custom brain parcellations by combining multiple atlases for different brain regions.{bcolors.ENDC}
"""

    epilog = f"""
{bcolors.BOLD}Examples:{bcolors.ENDC}

  {bcolors.OKYELLOW}# List all available parcellations per supra-region{bcolors.ENDC}
  chimera --regions

  {bcolors.OKYELLOW}# Basic usage with a 10-character parcellation code{bcolors.ENDC}
  chimera -b /data/study -d /data/study/derivatives -p HFMIIIFIFN

  {bcolors.OKYELLOW}# Interactive code generator (builds --parcodes, --seg, --scale, --growwm){bcolors.ENDC}
  chimera -b /data/study -d /data/study/derivatives -g

  {bcolors.OKYELLOW}# Multi-scale cortical parcellation (Schaefer 400 parcels, 7 networks){bcolors.ENDC}
  chimera -b /data/study -d /data/study/derivatives -p SFMIIIFIF -s 400 -e 7n

  {bcolors.OKYELLOW}# Grow GM labels 0 and 2 mm into white matter{bcolors.ENDC}
  chimera -b /data/study -d /data/study/derivatives -p HFMIIIFIFN -gw 0,2

  {bcolors.OKYELLOW}# Process specific subjects using 8 parallel threads{bcolors.ENDC}
  chimera -b /data/study -d /data/study/derivatives -p HFMIIIFIFN \\
          -ids sub-001,sub-002,sub-003 -n 8

  {bcolors.OKYELLOW}# Specify a custom FreeSurfer subjects directory{bcolors.ENDC}
  chimera -b /data/study -d /data/study/derivatives -p HFMIIIFIFN \\
          -fr /data/study/derivatives/freesurfer

{bcolors.BOLD}For more information:{bcolors.ENDC}
  Use --regions to see all available parcellation codes for each brain region.
"""

    p = argparse.ArgumentParser(
        formatter_class=ColoredHelpFormatter, description=description, epilog=epilog
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
        "-gw",
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
        "--gencode",
        "-g",
        action="store_true",
        required=False,
        help=f"{bcolors.BOLD}Launch interactive code generator{bcolors.ENDC}\n"
        f"Start the CHIMERA parcellation code generator to build --parcodes\n"
        f"(and optionally --seg / --scale) interactively.\n",
        default=False,
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

    # Show help when called with no arguments at all
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(0)

    # Launch the interactive code generator to fill --parcodes / --seg / --scale / --growwm
    if args.gencode:
        from chimera.chimera_code_generator import (
            _prompt_growwm as _cg_growwm,
        )
        from chimera.chimera_code_generator import (
            build_code_string as _cg_code_str,
        )
        from chimera.chimera_code_generator import (
            load_parcellations_info as _cg_load,
        )
        from chimera.chimera_code_generator import (
            run_interactive as _cg_run,
        )

        _cg_data = _cg_load()
        try:
            _cg_result = _cg_run(_cg_data)
            _cg_growwm_vals = _cg_growwm()
        except (KeyboardInterrupt, EOFError):
            print("\n  Interrupted -- no code generated.")
            sys.exit(0)

        _cg_code_string = _cg_code_str(_cg_result)
        args.parcodes = [_cg_code_string]

        # Collect seg and scale values from selected parcels
        import re as _re

        _all_segs, _all_scales = [], []
        for _choice in _cg_result.values():
            for _p in _choice.get("selected_parcels") or []:
                _m_seg = _re.search(r"seg-([^_]+)", _p)
                _m_scale = _re.search(r"scale-([^_]+)", _p)
                if _m_seg and _m_seg.group(1) not in _all_segs:
                    _all_segs.append(_m_seg.group(1))
                if _m_scale and _m_scale.group(1) not in _all_scales:
                    _all_scales.append(_m_scale.group(1))
        if _all_segs:
            args.seg = [",".join(_all_segs)]
        if _all_scales:
            args.scale = [",".join(_all_scales)]
        if _cg_growwm_vals != ["0"]:
            args.growwm = [",".join(_cg_growwm_vals)]

    global bids_dirs, supra_dict, deriv_dirs, fssubj_dirs, parcodes, pipe_json

    pipe_json = args.config

    if isinstance(args.config, list):
        pipe_json = args.config[0]

    if args.regions is True:
        print("\n")
        mess = "Available parcellations for each supra-region"
        print(
            f"{bcolors.BOLD}{bcolors.PURPLE}{mess}{bcolors.ENDC}{bcolors.ENDC}: "
        )
        _print_availab_parcels()
        sys.exit()

    # --bidsdir and --parcodes are required for an actual run
    if args.bidsdir is None or args.parcodes is None:
        print("--bidsdir and --parcodes are REQUIRED arguments")
        p.print_help()
        sys.exit()

    bids_dirs = args.bidsdir[0].split(sep=",")
    # Remove possible empty elements
    bids_dirs = [x for x in bids_dirs if x]

    for bids_dir in bids_dirs:
        if not os.path.isdir(bids_dir):
            print("Please, supply a valid BIDs directory.")
            print(
                f"The supplied BIDs directory does not exist: {bcolors.BOLD}{bcolors.OKRED}{bids_dir}{bcolors.ENDC}{bcolors.ENDC}: is not supplied. "
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
    """A simple progress indicator for the concurrent futures.

    :param future: future object
    """
    global lock, n_subj, n_comp, pb, pb1, chim_code
    # obtain the lock
    with lock:
        # update the counter
        n_comp += 1
        # report progress
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
    growwm: list = None,
    mixwm: bool = False,
    nthreads: int = 1,
):
    """Prepare chimera to build the parcellations.

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

    # Use a fresh list rather than a mutable default argument
    if growwm is None:
        growwm = ["0"]

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

            # Provide the resolved pipeline configuration to the object (replaces
            # the module-level global the original implementation relied on).
            chim_obj.pipe_dict = pipe_dict

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

            if nthreads == 1:

                pb1 = pb.add_task(
                    f"[red]Processing: Subject ({1}/{n_subj}) ", total=n_subj
                )
                for i, t1 in enumerate(t1s):

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
                    # send in the tasks, keeping the mapping to the subject so
                    # failures can be reported afterwards.
                    future2t1 = {
                        executor.submit(
                            chim_obj.build_parcellation,
                            t1s[i],
                            bids_dir,
                            deriv_dir,
                            fssubj_dir,
                            growwm,
                            mixwm,
                        ): t1s[i]
                        for i in range(n_subj)
                    }

                    # register the progress indicator callback
                    for future in future2t1:
                        future.add_done_callback(progress_indicator)

                # The ThreadPoolExecutor context manager waits for all tasks to
                # complete. Collect and report any subjects that failed instead of
                # silently swallowing the exception.
                for future, t1 in future2t1.items():
                    exc = future.exception()
                    if exc is not None:
                        failed.append((t1, exc))

                if failed:
                    print(
                        f"\n[{chim_code}] {len(failed)} subject(s) failed to process:"
                    )
                    for t1, exc in failed:
                        print(f"  - {os.path.basename(t1)}: {exc}")

                elapsed_time = time.perf_counter() - start_time
                print(
                    f"[{chim_code}] Processed {n_subj} subject(s) in "
                    f"{elapsed_time:.1f} s."
                )


def main():
    # 0. Handle inputs
    parser = _build_args_parser()
    args = parser.parse_args()

    print(args)
    if args.verbose is not None:
        v = int(args.verbose[0])
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
