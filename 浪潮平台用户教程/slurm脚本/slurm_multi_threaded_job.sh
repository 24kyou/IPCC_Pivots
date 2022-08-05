#!/bin/bash
#SBATCH --job-name=parallel_job              # Job name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yangyanwei@inspur.com    # Where to send mail	
#SBATCH --nodes=1                            # Run all processes on a single node	
#SBATCH --ntasks=1                           # Run a single task		
#SBATCH --cpus-per-task=4                    # Number of CPU cores per task
#SBATCH --mem=1gb                            # Job memory request
#SBATCH --time=00:05:00                      # Time limit hrs:min:sec
#SBATCH --output=parallel_%j.log             # Standard output and error log
pwd; hostname; date

echo "Running application program on $SLURM_CPUS_ON_NODE CPU cores"

/data/training/SLURM/parallel_job/application

date
