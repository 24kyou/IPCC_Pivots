#!/bin/bash
#SBATCH --nodelist=n00003
#SBATCH --time=10:00:10
#SBATCH --ntasks-per-node=56
#SBATCH --partition=G1Part
#SBATCH -J vasp                      # 作业名为 test

echo $SLURM_NPROCS

source  /es01/software/intel2020/compilers_and_libraries_2020.4.304/linux/bin/compilervars.sh intel64

mpirun -np $SLURM_NPROCS ./vasp_std

