#!/bin/bash -l

# Scheduling PARAMETERS

# Load R module:
## module load R 

# Name of the job - you'll probably want to customize this.
#SBATCH --job-name=poc-03
# Tell Gauss how much memory per CPU your job will use:
#SBATCH --mem=22120

# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o /home/cconley/scratch-data/neta-poc/mrna-tuning/03/log.out
#SBATCH -e /home/cconley/scratch-data/neta-poc/mrna-tuning/03/log.err

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=31

#Project Directory
#SBATCH -D /home/cconley/repos/neta-poc/mrna-tuning/

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
/usr/bin/R --no-save --no-restore --no-site-file --no-init-file --args 31 < /home/cconley/repos/neta-poc/mrna-tuning/03-cvvote-spacemap.R
