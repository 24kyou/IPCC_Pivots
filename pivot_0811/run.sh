#!/bin/bash
#SBATCH --job-name=Tt0811
#SBATCH --partition=IPCC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --error=%j.err
#SBATCH --output=%j.out


source /public1/soft/oneAPI/2022.1/setvars.sh

export OMP_NUM_THREADS=32
export OMP_PROC_BIND=true
export OMP_PLACES=cores

icc pivot_0811.c -Ofast -ipo -qopenmp -ffast-math  -o pivot_0811
./pivot_0811
