#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --time=100:00:00
#SBATCH --partition=D1Part_test
##SBATCH --nodelist=n08321
#SBATCH -J abinit
#SBATCH -o abinit.%j.out
srun hostname >./hostfile
echo $SLURM_NTASKS
source /es01/software/profile.d/intel2020.sh
inputfile=btm.files
#mpirun -genv I_MPI_FABRICS=shm:ofi -machinefile ./hostname -np $SLURM_NTASKS /es01/software/apps/abinit-8.10.3/bin/abinit  <$inputfile |tee log.$$
mpirun -genv I_MPI_FABRICS=shm:ofi -machinefile hostfile -np $SLURM_NTASKS /es01/software/apps/abinit-8.10.3/bin/abinit  <$inputfile
