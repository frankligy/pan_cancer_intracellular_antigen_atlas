#!/bin/bash
#SBATCH --partition=gpu4_medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=150Gb
#SBATCH --job-name="16s"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:v100:1


# module load sratoolkit/2.10.7


# declare -a samples=(
#     SRR24707480
#     SRR24707566
#     SRR24707581
#     SRR24708056
#     SRR24708135
#     SRR24708182
#     SRR24708286
#     SRR24707665
#     SRR24707772
#     SRR24708373
#     SRR24707522
#     SRR24707634
#     SRR24707753
#     SRR24708019
#     SRR24708153
#     SRR24708277
# )

# for f in ${samples[@]}; do echo $f; done | xargs -P 10 -I {} sh -c "fasterq-dump -e 20 {}" 
# # cut -d ',' -f 1 SraRunTable.csv | xargs -P 10 -I {} sh -c "fasterq-dump -e 20 {}" 
# for file in *.fastq; do echo $file; done | xargs -P 10 -I {} sh -c "gzip {}" 


# cut primer
# function remove_primer() {
#     SAMPLE=$1
#     ${CUTADAPT} --discard-untrimmed --pair-filter any -n 1 \
#                 -g ACTCCTACGGGAGGCAGCAG \
#                 -G GGACTACHVGGGTWTCTAAT \
#                 -o ${INTDIR}/trimmed/${SAMPLE}_trimmed_1.fastq \
#                 -p ${INTDIR}/trimmed/${SAMPLE}_trimmed_2.fastq \
#                 ${INTDIR}/${SAMPLE}_1.fastq ${INTDIR}/${SAMPLE}_2.fastq
# }
# module load parallel
# export -f remove_primer
# export TMPDIR=/gpfs/data/yarmarkovichlab/Frank/pan_cancer/pathogen_seq
# export CUTADAPT=/gpfs/data/yarmarkovichlab/russell_data/cutadapt_env/bin/cutadapt
# export INTDIR=/gpfs/data/yarmarkovichlab/Frank/pan_cancer/pathogen_seq/CX
# cat ${INTDIR}/samples.txt | parallel -P 1 remove_primer {} > ${INTDIR}/cutadapt_log.txt


# module load singularity/3.9.8
# https://github.com/blekhmanlab/docker_dada2/tree/main
# https://benjjneb.github.io/dada2/tutorial.html
# singularity build --sandbox ./DATA2/ docker://blekhmanlab/dada2 (but I think sif works as well, for sandbox, you need to mkdir /gpfs)
# singularity shell -B ./PORT:/mnt --writable ./DATA2/
# launch R

module load singularity/3.9.8
singularity run -B ./CX:/mnt --writable ./DATA2/ Rscript /mnt/rRNA.R

