#!/bin/sh
#SBATCH --job-name=mlgp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -p regular
#SBATCH --array=<START>-<END>%8
#SBATCH --time=03:30:00
index=$(printf "%01d" ${SLURM_ARRAY_TASK_ID})
module load intel_compilers
export OMP_NUM_THREADS=24
export CRDATA="/home/nicusor/software/cretin.v2_20/data"
/home/nicusor/software/cretin.v2_20_21_03_09/cretin.ifort run_${index}.gen > run_${index}.out
exit 0
