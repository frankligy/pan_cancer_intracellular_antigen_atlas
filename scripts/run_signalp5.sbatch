#!/bin/bash
#SBATCH --partition=a100_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=150Gb
#SBATCH --job-name="neoverse-pan"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:a100:1


export PATH=/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/signalp-5.0b/bin:$PATH

signalp -fasta /gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/summary/membrane_protein.fasta -org euk -format short -prefix run