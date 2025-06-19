#!/bin/bash
#SBATCH --partition=gpu4_medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=150Gb
#SBATCH --job-name="neoverse-pan"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:v100:1




function run_unmapped_kraken2 () {
    # download the kraken2 pre-built database: https://benlangmead.github.io/aws-indexes/k2 using PlusPF 
    # curl -o k2_pluspf_20230605.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20230605.tar.gz
    # and only need opts.k2d, hash.k2d, taxo.k2d, name it as kraken2_db

    if [ ! -d ${PATHOGEN_TMP} ]; then
        mkdir ${PATHOGEN_TMP}
    fi

    module load samtools
    samtools view -bh -f 12 ${STAR_TMP}/${SAMPLE}.rna_seq.genomic.gdc_realn.bam > ${PATHOGEN_TMP}/both_unmapped.bam

    module load samtools
    module load bedtools

    BAM=${PATHOGEN_TMP}/both_unmapped.bam
    samtools sort -n $BAM -o ${PATHOGEN_TMP}/both_unmapped.qsort.bam
    bedtools bamtofastq -i ${PATHOGEN_TMP}/both_unmapped.qsort.bam -fq ${PATHOGEN_TMP}/both_unmapped_1.fastq -fq2 ${PATHOGEN_TMP}/both_unmapped_2.fastq


    module load kraken/2.0.8
    kraken2 --db ${KRAKEN2_DB} \
            --threads 12 \
            --paired \
            --classified-out ${PATHOGEN_TMP}/test_cseqs#.fq \
            --unclassified-out ${PATHOGEN_TMP}/test_ucseqs#.fq \
            --report ${PATHOGEN_TMP}/test_report.txt \
            --use-names \
            --use-mpa-style \
            --output ${PATHOGEN_TMP}/test_output.txt \
            --confidence 0.9 \
            ${PATHOGEN_TMP}/both_unmapped_1.fastq ${PATHOGEN_TMP}/both_unmapped_2.fastq
    module purge
    
}

KRAKEN2_DB=/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/pathogen/kraken2_db
STAR_TMP=/gpfs/data/yarmarkovichlab/Frank/pan_cancer/n_circulans_investigate

while read line;do
    arg1=$(echo "${line}" | awk -F "\t" '{print $11}')
    arg2=$(echo "${line}" | awk -F "\t" '{print $3}')
    if [[ "$arg1" != "sample_id" ]]; then
        PATHOGEN_TMP=/gpfs/data/yarmarkovichlab/Frank/pan_cancer/n_circulans_investigate/${arg1}
        arr_arg2=(${arg2//./ })
        SAMPLE=${arr_arg2[0]}
        echo ${PATHOGEN_TMP}
        echo ${SAMPLE}
        run_unmapped_kraken2
    fi
done < manifest.txt