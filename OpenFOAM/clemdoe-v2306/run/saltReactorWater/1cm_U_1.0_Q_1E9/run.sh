#!/bin/bash
# File: run.sh
# Author: Cristian Garrido Tamm (cristian.garrido@idom.com)
# Purpose:
# Creation Date: Mon 15 Jan 2024 05:42:54 PM CET
# Last Modified: Sat 09 Mar 2024 09:14:29 PM CET
######################
### SBATCH OPTIONS ###
######################

# Job Name
#SBATCH -J 33_100cm_U_0.25_Q_100E6 

# Name of stdout file
#SBATCH --output="run-%j.out"

# Name of stderr file
#SBATCH --error="run-%j.err"

# Calculation Time d-hh:mm:ss
# #SBATCH --time=7-00:00:00 

# Number of Nodes
#SBATCH --nodes=1

# Number of tasks
#SBATCH --ntasks=40

# Number of CPUs for each task
#SBATCH --cpus-per-task=1

# List of nodes to use
# #SBATCH --nodelist=[c01, c02, c03]

# How much memory you need.
# --mem will define memory per node and
# --mem-per-cpu will define memory per CPU/core. Choose one of those.
# #SBATCH --mem-per-cpu=7900

# Send e-mail at the end (Not working!)
#SBATCH --mail-type=end
#SBATCH --mail-user=cristian.garrido@idom.com

set -eu

# . /home/shared/spack/share/spack/setup-env.sh

spack load openfoam@2306

source ~/.venvs/pyfoam/bin/activate

./Allrun-parallel

pyFoamPlotWatcher.py --with-all --implementation='matplotlib' --hardcopy --solver-not-running-anymore --progress log.chtMultiRegionSimpleFoam

pyFoamClearCase.py --no-allclean-script --keep-last --keep-postprocessing --processors-remove --vtk-keep .

