import os
import subprocess
import uuid
from glob import glob
from pathlib import Path

import clabtoolkit.misctools as cltmisc


def run_command(
    bash_args: list, cont_tech: str = "local", cont_image: str = None
) -> str:
    """
    Run a command locally or inside a container and return its standard output.

    The command runs without a shell, so arguments are passed literally: globs
    and environment variables are not expanded. Wrap the command in
    ``["bash", "-c", "..."]`` when shell expansion is needed.

    Parameters:
    ----------
    bash_args : list
        Command and its arguments.

    cont_tech : str
        Container technology (e.g. singularity, docker or local).

    cont_image : str
        Container image.

    Returns:
    --------
    stdout : str
        Standard output of the command.

    Raises:
    -------
    RuntimeError
        If the command exits with a non-zero status. The message includes the
        command and its standard error.
    """
    cmd_cont = cltmisc.generate_container_command(bash_args, cont_tech, cont_image)
    result = subprocess.run(cmd_cont, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit status {result.returncode}: "
            f"{' '.join(cmd_cont)}\n{result.stderr.strip()}"
        )
    return result.stdout


def launch_fsl_first(
    t1: str,
    first_parc: str,
    cont_tech: str = "local",
    cont_image: str = None,
    overwrite=False,
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

    overwrite : bool
        Overwrite the results.

    Returns:
    --------

    Raises:
    -------
    RuntimeError
        If FIRST, or moving its output into place, fails.
    """
    fsl_outdir = os.path.dirname(first_parc)

    if not os.path.isfile(first_parc) or overwrite:
        Path(fsl_outdir).mkdir(parents=True, exist_ok=True)

        # A unique prefix keeps runs that share an output directory (e.g. two
        # T1w images of the same subject processed in parallel) from reading or
        # deleting each other's intermediate files.
        tmp_prefix = os.path.join(fsl_outdir, f"tmp_first_{uuid.uuid4().hex}")

        try:
            run_command(
                ["run_first_all", "-i", t1, "-o", tmp_prefix], cont_tech, cont_image
            )
            run_command(
                ["mv", tmp_prefix + "_all_fast_firstseg.nii.gz", first_parc],
                cont_tech,
                cont_image,
            )
        finally:
            # FIRST writes <prefix>_* files and a <prefix>.logs directory. The
            # paths are listed here, rather than globbed by the command, because
            # commands run without a shell. The removal runs through the same
            # container as FIRST, so files the container created are removable.
            tmp_files = glob(tmp_prefix + "*")
            if tmp_files:
                run_command(["rm", "-rf"] + tmp_files, cont_tech, cont_image)
