#!/bin/bash
#SBATCH --job-name=parallel_job_test         # Job name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yangyanwei@inspur.com    # Where to send mail	
#SBATCH --nodes=1                            # Run all processes on a single node	
#SBATCH --ntasks=4                           # Number of processes
#SBATCH --mem=1gb                            # Total memory limit
#SBATCH --time=01:00:00                      # Time limit hrs:min:sec
#SBATCH --output=multiprocess_%j.log         # Standard output and error log
date;hostname;pwd

module load python/3

python your_script.py

date
