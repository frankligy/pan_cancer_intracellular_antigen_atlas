#!/bin/bash
#SBATCH --partition=a100_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-00:00:00
#SBATCH --mem=50Gb
#SBATCH --job-name="data_transfer"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:a100:1



# PXD027766 just plasma, don't do pleural effusion or cells
declare -a samples=(
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/05/PXD027766/Seq62941_QE3.raw
)

for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P ../immuno/dummy {}"


