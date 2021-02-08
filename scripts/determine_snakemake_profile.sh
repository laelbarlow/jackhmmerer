#!/usr/bin/env bash

# Determine what string to use to define Snakemake profile in current
# environment.

snakemake_profile=""

if test "$(command -v sbatch)"; then
    echo SLURM job scheduler detected.
    snakemake_profile="slurm_jackhmmerer"

elif qstat -help; then
    printf "(SUN/Univa) Grid Engine job scheduler detected."
    snakemake_profile="sge_jackhmmerer"

elif pbs-config --version; then
    printf "PBS-TORQUE job scheduler detected."
    snakemake_profile="pbs-torque_jackhmmerer"

else
    echo Unable to detect job scheduler.
    #exit 1

fi
