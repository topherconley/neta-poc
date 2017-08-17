#!/bin/bash -l

# Scheduling PARAMETERS

# Load R module:
## module load R 

# Name of the job - you'll probably want to customize this.
#SBATCH --job-name=poc-02
# Tell Gauss how much memory per CPU your job will use:
#SBATCH --mem=22120

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o /home/cconley/scratch-data/neta-poc/tuning/02/log.out
#SBATCH -e /home/cconley/scratch-data/neta-poc/tuning/02/log.err

#SBATCH --nodes=8
#SBATCH --ntasks-per-node=15

#Project Directory
#SBATCH -D /home/cconley/repos/neta-poc/tuning/

###############################################################################
# Bash options
# End the job at the first sign of an error
set -e
#Print each command to stdout before executing it
set -v
#tell me the host of my job
hostname

# Execute each of the jobs with a different index (the R script will then process
# this to do something different for each index):
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file --args 15 < /home/cconley/repos/neta-poc/tuning/02-cvvote-space.R
