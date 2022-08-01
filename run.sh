#!/bin/bash
#SBATCH --job-name=pivot0801
#SBATCH --partition=IPCC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=%j.err
#SBATCH --output=%j.out

gcc pivot0801.c -lm -O2 -fopenmp -o pivotest
./pivotest

