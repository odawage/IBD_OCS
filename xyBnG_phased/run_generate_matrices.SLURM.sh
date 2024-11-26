#!/bin/bash
#SBATCH --job-name=gen_mat     # Job name
#SBATCH --output=gen_mat_%j.out   # Output file
#SBATCH --nodes=1
#SBATCH --ntasks=1                  # Number of tasks (cores)
#SBATCH --cpus-per-task=29
#SBATCH --partition=hugemem-avx2,orion,hugemem   # Partition name
#SBATCH --mem=1400G                     # 400 Memory required per node
#SBATCH --exclude=cn-14,cn-11

module load Julia/1.10.0-linux-x86_64
#module load Python/3.7.4-GCCcore-8.3.0  
#module load Miniconda3

module list 

eval "$(conda shell.bash hook)"

##Load conda environment 
conda activate simulation

# # export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
# cd /mnt/SCRATCH/odwa/xyBnG_sim

# cp /mnt/users/odwa/paper-2/xyBnG_phased/running_sim.jl .

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --threads ${SLURM_CPUS_PER_TASK} generate_matrices.jl

#cp rst/cattle $HOME/paper-2/xyBnG_phased