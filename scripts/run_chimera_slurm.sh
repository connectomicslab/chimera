#!/usr/bin/env bash
#
# Run the Chimera parcellation tool as a SLURM array job.
#
# Usage:
#   sbatch --array=1-<n> run_chimera_slurm.sh BIDS_DIR DERIV_DIR IDS_FILE [extra chimera args...]
#
#   BIDS_DIR   BIDS dataset directory
#   DERIV_DIR  derivatives directory
#   IDS_FILE   text file with one ID per line; line N is processed by array task N.
#              An ID is the T1w file name without suffix/extension, e.g.
#              sub-137S4816_ses-m03_acq-3T_run-01
#
# Everything else is configured through environment variables (see the block
# below) or by appending flags after IDS_FILE, which are forwarded verbatim to
# chimera and win over the environment variables.
#
# Examples:
#   # defaults
#   sbatch --array=1-50 run_chimera_slurm.sh /data/bids /data/deriv ids.txt
#
#   # different code + scale for one study, no edits to this file
#   PARCODES=HFMIIIFIFN SCALE=3 sbatch --array=1-50 run_chimera_slurm.sh \
#       /data/bids /data/deriv ids.txt
#
#   # anything chimera accepts, passed straight through
#   sbatch --array=1-50 run_chimera_slurm.sh /data/bids /data/deriv ids.txt \
#       --seg kong7n --mergectx
#
#   # run a single subject outside SLURM
#   TASK_ID=3 ./run_chimera_slurm.sh /data/bids /data/deriv ids.txt

set -euo pipefail

#########################################################################################
########################################## Settings #####################################
#########################################################################################

# Every value below can be overridden from the environment:
#   PARCODES=... GROWWM=... sbatch ... run_chimera_slurm.sh ...
# Leave one empty to omit the corresponding chimera flag entirely.

chim_config="${CHIMERA_CONFIG-/cluster/home/yalemang/Worktools/cluster_scripts/pipe_config.json}"
parcodes="${PARCODES-LFIIHISIFN}"
growwm="${GROWWM-1}"
scale="${SCALE:-}"
seg="${SEG:-}"
nthreads="${NTHREADS-1}"
verbose="${VERBOSE:-}"
mergectx="${MERGECTX:-0}"   # 1 adds --mergectx
overwrite="${FORCE:-0}"         # 1 adds --overwrite

# FreeSurfer SUBJECTS_DIR. Defaults to <DERIV_DIR>/freesurfer.
fs_dir="${FREESURFER_DIR:-}"

# Conda environment holding chimera. Set CONDA_SH='' to skip activation
# (e.g. when the environment is already active or chimera is on the PATH).
conda_sh="${CONDA_SH-/data/hagmann_group/install/miniconda3/etc/profile.d/conda.sh}"
conda_env="${CONDA_ENV:-/data/hagmann_group/conda_environments/chimera-env}"

# Name of this pipeline step, used for the lock and error directories.
step_name="${STEP_NAME:-step6-chimera}"

# Set SKIP_LOCK=1 to ignore the lock files and always (re)process.
skip_lock="${SKIP_LOCK:-0}"

#########################################################################################
########################################## Main program #################################
#########################################################################################

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 BIDS_DIR DERIV_DIR IDS_FILE [extra chimera args...]" >&2
    exit 2
fi

bids_dir=$1
out_dir=$2
ids_file=$3
shift 3
extra_args=("$@")   # forwarded verbatim to chimera, overrides the settings above

for path in "$bids_dir" "$ids_file"; do
    if [ ! -e "$path" ]; then
        echo "ERROR: no such file or directory: $path" >&2
        exit 2
    fi
done

# Array task index; TASK_ID lets the script run outside of SLURM.
task_id="${SLURM_ARRAY_TASK_ID:-${TASK_ID:-}}"
if [ -z "$task_id" ]; then
    echo "ERROR: neither SLURM_ARRAY_TASK_ID nor TASK_ID is set." >&2
    echo "       Submit with 'sbatch --array=1-<n>' or set TASK_ID=<line number>." >&2
    exit 2
fi

subj_full_id=$(sed -n "${task_id}p" "$ids_file")
if [ -z "$subj_full_id" ]; then
    echo "ERROR: line $task_id of $ids_file is empty or out of range." >&2
    exit 2
fi

[ -n "$fs_dir" ] || fs_dir="$out_dir/freesurfer"

lock_dir="$out_dir/processing-lockfiles/$step_name"
error_dir="$out_dir/pipeline-errors/$step_name"
mkdir -p "$lock_dir" "$error_dir"

lockfile="$lock_dir/${subj_full_id}_desc-${step_name//-/}_LOCK.lock"
if [ "$skip_lock" != "1" ] && [ -f "$lockfile" ]; then
    echo "Already processed (lock file present): $lockfile"
    exit 0
fi

#  ------- Build the chimera command
# Flags are only added when the corresponding setting is non-empty, so an
# unwanted option is switched off with e.g. GROWWM=''.
chimera_args=(
    --bidsdir "$bids_dir"
    --derivdir "$out_dir"
    --freesurferdir "$fs_dir"
    -ids "$subj_full_id"
    --parcodes "$parcodes"
    --nthreads "$nthreads"
)
if [ -n "$chim_config" ]; then chimera_args+=(--config "$chim_config"); fi
if [ -n "$growwm" ];      then chimera_args+=(--growwm "$growwm");      fi
if [ -n "$scale" ];       then chimera_args+=(--scale "$scale");        fi
if [ -n "$seg" ];         then chimera_args+=(--seg "$seg");            fi
if [ -n "$verbose" ];     then chimera_args+=(--verbose "$verbose");    fi
if [ "$mergectx" = "1" ]; then chimera_args+=(--mergectx);              fi
if [ "$overwrite" = "1" ];    then chimera_args+=(--overwrite);                 fi

# The '+' expansion keeps an empty array from tripping 'set -u' on bash < 4.4.
chimera_args+=(${extra_args[@]+"${extra_args[@]}"})

echo "Processing        : $subj_full_id  (task $task_id)"
echo "BIDS_DIR          : $bids_dir"
echo "DERIVATIVES_DIR   : $out_dir"
echo "FREESURFER_DIR    : $fs_dir"
echo "CONFIG_FILE       : ${chim_config:-<none>}"
echo " "
echo "------ Command -----------"
echo " "
echo " 1. Running chimera parcellation"
echo "------------------------------------"
printf 'chimera'; printf ' %q' "${chimera_args[@]}"; printf '\n'
echo " "

#  ------- Loading modules
# module load singularityce
# module load miniconda3

if [ -n "$conda_sh" ]; then
    # conda's own scripts are not written for 'set -u'.
    set +u
    # shellcheck disable=SC1090
    source "$conda_sh"
    conda activate "$conda_env"
    set -u
fi

# Take the lock now, so a concurrent submission does not process the same
# subject, and release it again if chimera fails so the job can be retried.
touch "$lockfile"

status=0
chimera "${chimera_args[@]}" || status=$?

if [ "$status" -ne 0 ]; then
    rm -f "$lockfile"
    echo "$subj_full_id exited with status $status on $(date -Is)" \
        >> "$error_dir/${subj_full_id}_desc-${step_name//-/}_ERROR.txt"
    echo "ERROR: chimera failed for $subj_full_id (status $status)" >&2
    exit "$status"
fi

echo "==================="
echo " "
