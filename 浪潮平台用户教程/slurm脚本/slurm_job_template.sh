#!/bin/bash
#SBATCH --job-name="testJob"
#SBATCH --partition=computerPartiton
#SBATCH --qos=qos_large
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --output=job.%j.out

echo "current job: ${SLURM_JOB_ID}:${SLURM_JOB_NAME}"
echo "  running on: ${SLURM_JOB_PARTITION}:${SLURM_JOB_NODELIST}"
echo "  allocated resources: CPU ${SLURM_CPUS_ON_NODE} / GPU ${SLURM_GPUS_PER_NODE}"
sleep 60

