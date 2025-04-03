#!/bin/bash 

#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=5TPN
#SBATCH --time=48:00:00
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A30:1
#SBATCH -o logs/5TPN.output
#SBATCH -e logs/5TPN.error


input_path=/vast/scratch/users/$USER/inputs
output_path=/vast/scratch/users/$USER/outputs
schedule_path=/vast/scratch/users/$USER/schedule

module load apptainer/1.3.3

##comment the next three files if already created.
mkdir $input_path
mkdir $output_path
mkdir $schedule_path

wget -P $input_path https://files.rcsb.org/view/5TPN.pdb

apptainer run --nv \
-B /vast/projects/alphafold/databases/predictedDb/RFDiffusion/models:$HOME/models  \
-B $input_path:$HOME/inputs \
-B $output_path:$HOME/outputs \
-B $schedule_path:$HOME/schedule \
/stornext/System/data/apps/rc-tools/rc-tools-1.0/bin/tools/rfdiffusion.sif \
inference.output_prefix=$HOME/outputs/ \
inference.model_directory_path=$HOME/models \
inference.input_pdb=$HOME/inputs/5TPN.pdb \
inference.schedule_directory_path=$HOME/schedule \
inference.num_designs=3 \
'contigmap.contigs=[10-40/A163-181/10-40]'
