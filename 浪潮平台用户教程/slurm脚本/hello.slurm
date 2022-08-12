#!/bin/bash
#SBATCH --job-name=mpi_job_test              # Job name
#SBATCH --output=testSlurmJob.%j.out    # Stdout (%j expands to jobId)
#SBATCH --error=testSlurmJob.%j.err     # Stderr (%j expands to jobId)
#SBATCH --ntasks=24                          # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1                    # Number of cores per MPI task 
#SBATCH --nodes=2                            # Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=12                 # Maximum number of tasks on each node
#SBATCH --ntasks-per-socket=6                # Maximum number of tasks on each socket
#SBATCH --mem-per-cpu=600mb                  # Memory (i.e. RAM) per processor
#SBATCH --time=00:05:00                      # Wall time limit (days-hrs:min:sec)
srun hostname >./hostfile
echo $SLURM_NTASKS
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo $SLURM_NPROCS
echo $SLURM_NPROCS
source /software/intel2018/compilers_and_libraries_2018/linux/bin/compilervars.sh intel64
source /software/intel2018/mkl/bin/mklvars.sh intel64
source /software/intel2018/impi/2018.1.163/bin64/mpivars.sh
mpirun -genv I_MPI_FABRICS=shm:ofa -machinefile hostfile -np $SLURM_NTASKS  /home/test/test/hello
sleep 10
