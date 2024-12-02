#!/bin/bash
#SBATCH --job-name=sim_cow     # Job name
#SBATCH --output=sim_cow_%j.out   # Output file
#SBATCH --nodes=1
#SBATCH --ntasks=1                  # Number of tasks (cores)
#SBATCH --cpus-per-task=29
#SBATCH --partition=hugemem-avx2,orion,hugemem   # Partition name
#SBATCH --mem=100G                     # 400 Memory required per node
#SBATCH --exclude=cn-14,cn-11

module load Julia/1.10.0-linux-x86_64
#module load Python/3.7.4-GCCcore-8.3.0  
#module load Miniconda3

module list 

eval "$(conda shell.bash hook)"

##Load conda environment 
conda activate simulation


###Start the job by checking if the $TMPDIR is available

## Copying data to local node for faster computation

cd $TMPDIR

#Check if $USER exists in $TMPDIR

if [[ -d $USER ]]
        then
                echo "$USER exists on $TMPDIR"
        else
                mkdir $USER
fi


echo "copying files to" $TMPDIR/$USER

##Create a tempdirectory specific for each job

cd $USER
mkdir tmpDir_of.$SLURM_JOB_ID
cd tmpDir_of.$SLURM_JOB_ID

echo "My working directory is:"
pwd

#Rsync the fastq files and the reference

##The following variable helps not to write all the rsync arguments each time

RSYNC='rsync -aL --no-perms --no-owner --no-group --no-t'

if [ -e /mnt/project/SimData/xyBnG/tskit_50/founder_files.tar.gz ]; then
    echo "The founder files are zipped and ready"
    mkdir -p tskit_50
    $RSYNC /mnt/project/SimData/xyBnG/tskit_50/founder_files.tar.gz ./tskit_50
    tar -xzvf ./tskit_50/founder_files.tar.gz -C ./tskit_50
        echo "files are zippd, moved and here"

else
    if [ -e /mnt/project/SimData/xyBnG/tskit_50/BosTau.lmp ]; then
        tar -cvf - -C /mnt/project/SimData/xyBnG/tskit_50/ BosTau.lmp BosTau.xy desc.txt | pigz --fast > /mnt/project/SimData/xyBnG/tskit_50/founder_files.tar.gz
        $RSYNC /mnt/project/SimData/xyBnG/tskit_50/founder_files.tar.gz ./tskit_50
        tar -xzvf ./tskit_50/founder_files.tar.gz -C ./tskit_50

        echo "files are zippd, moved and here"

    else
        echo "No founder simulation done, rerunning it"
    fi
fi


export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

julia --threads ${SLURM_CPUS_PER_TASK}  /mnt/users/odwa/IBD_OCS/xyBnG_phased/running_sim.jl


#cp rst/cattle $HOME/paper-2/xyBnG_phased

cd $TMPDIR/$USER/

$RSYNC tmpDir_of.$SLURM_JOB_ID /mnt/project/SimData/xyBnG/ 

rm -r tmpDir_of.$SLURM_JOB_ID 