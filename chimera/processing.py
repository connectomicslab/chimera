import os
import subprocess
from pathlib import Path


def launch_fsl_first(
    t1: str,
    first_parc: str,
    cont_tech: str = "local",
    cont_image: str = None,
    force=False,
):
    """
    This function executes the FIRST subcortical parcellation.

    Parameters:
    ----------
    t1 : str
        T1-weighted image filename.

    first_parc : str
        Ouput name for the resulting parcellation.

    cont_tech : str
        Container technology (e.g. singularity, docker or local).

    cont_image : str
        Container image.

    force : bool
        Overwrite the results.

    Returns:
    --------

    """
    fsl_outdir = os.path.dirname(first_parc)

    if not os.path.isfile(first_parc) or force:
        fsl_outdir = Path(fsl_outdir)
        fsl_outdir.mkdir(parents=True, exist_ok=True)

        cmd_bashargs = [
            "run_first_all",
            "-i",
            t1,
            "-o",
            str(fsl_outdir) + os.path.sep + "temp",
        ]
        cmd_cont = cltmisc.generate_container_command(
            cmd_bashargs, cont_tech, cont_image
        )  # Generating container command
        subprocess.run(
            cmd_cont, stdout=subprocess.PIPE, universal_newlines=True
        )  # Running container command

        cmd_bashargs = [
            "mv",
            os.path.join(str(fsl_outdir), "temp_all_fast_firstseg.nii.gz"),
            first_parc,
        ]
        cmd_cont = cltmisc.generate_container_command(
            cmd_bashargs, cont_tech, cont_image
        )  # Generating container command
        subprocess.run(
            cmd_cont, stdout=subprocess.PIPE, universal_newlines=True
        )  # Running container command

        cmd_bashargs = ["rm", "-rf", "temp*"]
        cmd_cont = cltmisc.generate_container_command(
            cmd_bashargs, cont_tech, cont_image
        )  # Generating container command
        subprocess.run(
            cmd_cont, stdout=subprocess.PIPE, universal_newlines=True
        )  # Running container command
