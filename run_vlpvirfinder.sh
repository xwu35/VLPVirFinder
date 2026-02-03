#!/bin/bash
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH -J run_vlpvirfinder
#SBATCH --output=slurm-output/slurm-%x.%j.out
#SBATCH --error=slurm-output/slurm-%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xiaofen

eval "$(conda shell.bash hook)"
conda activate snakemake

python VLPVirFinder.py \
    --reads_dir /path/to/dir/containing/sequences \
    --sample_info /path/to/sample_info_table \
    --output_dir output_dir 

