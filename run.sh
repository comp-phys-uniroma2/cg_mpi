#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=out.%j
#SBATCH --error=err.%j
#SBATCH --time=24:00:00
#SBATCH --partition=idra

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
module load gcc/gcc-11.1
module load mpi/openmpi-4.1.4-gcc11.1

srun  ./poisson_p 

