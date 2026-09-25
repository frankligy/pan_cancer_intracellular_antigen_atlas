#!/bin/bash
#SBATCH --partition=a100_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=1-00:00:00
#SBATCH --mem=150Gb
#SBATCH --job-name="neoverse-pan"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:a100:1


# /gpfs/data/yarmarkovichlab/softwares/NeoVerse_pan/integrate_tcga.py --config template_2.json

# /gpfs/data/yarmarkovichlab/softwares/NeoVerse_pan/NeoVerse_pan.py --config template_2.json --running_mode non_interactive

# /gpfs/data/yarmarkovichlab/softwares/NeoVerse_pan/run_maxquant.sh template_2.json

/gpfs/data/yarmarkovichlab/softwares/NeoVerse_pan/post_maxquant.py --config template_2.json

# /gpfs/data/yarmarkovichlab/softwares/NeoVerse_pan/rederive.py --config template_2.json

# module load singularity/3.1
# /gpfs/data/yarmarkovichlab/softwares/NeoVerse_pan/launch_portal.py --config template_2.json --running_mode set_up_rescue