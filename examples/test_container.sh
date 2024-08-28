#!/bin/bash 

#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=5TPN
#SBATCH --time=48:00:00
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A30:1
#SBATCH -o logs/5TPN.output
#SBATCH -e logs/5TPN.error

module load apptainer/1.3.3
mkdir /vast/scratch/users/$USER/inputs
mkdir /vast/scratch/users/$USER/outputs
mkdir /vast/scratch/users/$USER/schedule
wget -P /vast/scratch/users/$USER/inputs https://files.rcsb.org/view/5TPN.pdb

apptainer run --nv \
-B /vast/projects/alphafold/databases/predictedDb/RFDiffusion/models:$HOME/models  \
-B /vast/scratch/users/$USER/inputs:$HOME/inputs \
-B /vast/scratch/users/$USER/outputs:$HOME/outputs \
-B /vast/scratch/users/$USER/schedule:$HOME/schedule \
../rfdiffusion.sif \
inference.output_prefix=outputs/ \
inference.model_directory_path=models \
inference.input_pdb=inputs/5TPN.pdb \
inference.schedule_directory_path=schedule \
inference.num_designs=3 \
'contigmap.contigs=[10-40/A163-181/10-40]'