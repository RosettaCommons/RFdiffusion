#!/bin/bash 
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --job-name=5TPN
#SBATCH --time=48:00:00
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A30:1
#SBATCH -o logs/CTR.output
#SBATCH -e logs/CTR.error


input_path=/vast/scratch/users/$USER/inputs
output_path=/vast/scratch/users/$USER/outputs
schedule_path=/vast/scratch/users/$USER/schedule

module load apptainer/1.3.3

apptainer run --nv \
-B /vast/projects/alphafold/databases/predictedDb/RFDiffusion/models:$HOME/models \
-B $input_path:$HOME/inputs \
-B $output_path:$HOME/outputs \
-B $schedule_path:$HOME/schedule \
/stornext/System/data/apps/rc-tools/rc-tools-1.0/bin/tools/rfdiffusion.sif \
inference.output_prefix=outputs/ \
inference.model_directory_path=models \
inference.input_pdb=inputs/6asy.pdb \
inference.schedule_directory_path=schedule \
inference.num_designs=1000 \
denoiser.noise_scale_ca=0 \
denoiser.noise_scale_frame=0 \
'contigmap.contigs=[A24-629/0 50-50]'
