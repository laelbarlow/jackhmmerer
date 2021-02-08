#!/bin/bash
# properties = {properties}

# Load modules.
module load  \
nixpkgs/16.09 \
intel/2018.3 \
openmpi/3.1.4 \

# Execute job.
{exec_job}
