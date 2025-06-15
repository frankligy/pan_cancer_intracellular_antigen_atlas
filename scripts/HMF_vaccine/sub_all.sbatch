#!/bin/bash
#SBATCH --partition=a100_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=1-00:00:00
#SBATCH --mem=150Gb
#SBATCH --job-name="neoverse"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:a100:1


# /gpfs/data/yarmarkovichlab/softwares/NeoVerse/preprocess.sh template.json

# step 1
# /gpfs/data/yarmarkovichlab/softwares/NeoVerse/NeoVerse.sh template.json

# step 2
# module load samtools
# module load htslib/1.3
# module load singularity/3.9.8
# export LD_LIBRARY_PATH=/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/lib:$LD_LIBRARY_PATH
# /gpfs/data/yarmarkovichlab/softwares/NeoVerse/NeoVerse.py --config template.json --running_mode non_interactive

# # step 3
# /gpfs/data/yarmarkovichlab/softwares/NeoVerse/configure_maxquant.py --config chordoma_pdx.json 


# step 4
# /gpfs/data/yarmarkovichlab/softwares/NeoVerse/run_maxquant.sh template.json


/gpfs/data/yarmarkovichlab/softwares/NeoVerse/post_maxquant.py --config template.json

# /gpfs/data/yarmarkovichlab/softwares/NeoVerse/rederive.py --config template.json

# module load singularity/3.1
# /gpfs/data/yarmarkovichlab/softwares/NeoVerse/launch_portal.py --config template.json --running_mode generate_figures

# ./personal_mutation.py