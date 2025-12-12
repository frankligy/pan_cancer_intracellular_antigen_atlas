#!/bin/bash



function run_riborf_part1 () {


    gunzip -c ${FASTQ} > ${OUTDIR}/${SAMPLE}.fastq

    # module load perl/5.38.2
    # perl ${RIBORF_CODE}/RibORF.2.0/removeAdapter.pl -f ${OUTDIR}/${SAMPLE}.fastq -a ${ADAPTER_SEQ} -o ${OUTDIR}/adapter.${SAMPLE}.fastq
    # module purge

    ${CUTADAPT} -a ${ADAPTER_SEQ} -u 3 -m 15 -o ${OUTDIR}/adapter.${SAMPLE}.fastq ${OUTDIR}/${SAMPLE}.fastq
    

    module load bowtie2/2.3.1
    bowtie2-build ${RRNA_FASTA} ${OUTDIR}/hg19.ribosome
    bowtie2 -x ${OUTDIR}/hg19.ribosome -U ${OUTDIR}/adapter.${SAMPLE}.fastq --un ${OUTDIR}/norrna.adapter.${SAMPLE}.fastq -S ${OUTDIR}/ribosome.adapter.${SAMPLE}.sam
    rm ${OUTDIR}/hg19.ribosome*
    module purge


    module load tophat/2.1.1
    module load python/cpu/2.7.15
    module load bowtie2/2.3.1

    python $(which tophat) --GTF ${HG19_GTF} --no-convert-bam -o ${OUTDIR}/outputDir ${BOWTIE_HG19_INDEX_DIR}/hg19genome.index ${OUTDIR}/norrna.adapter.${SAMPLE}.fastq
    module purge


    module load perl/5.38.2
    module load r/4.3.2

    perl ${RIBORF_CODE}/RibORF.2.0/readDist.pl -f ${OUTDIR}/outputDir/accepted_hits.sam -g ${HG19_GENEPRED} -o ${OUTDIR}/outputDir 
    module purge

    rm ${OUTDIR}/${SAMPLE}.fastq
    rm ${OUTDIR}/adapter.${SAMPLE}.fastq
    rm ${OUTDIR}/norrna.adapter.${SAMPLE}.fastq
    rm ${OUTDIR}/ribosome.adapter.${SAMPLE}.sam


}


function run_riborf_part2 () {
    OFFSET_FILE_PATH=${OUTDIR}/outputDir/offset.correction.parameters.txt

    module load perl/5.38.2
    module load r/4.3.2
    perl ${RIBORF_CODE}/RibORF.2.0/offsetCorrect.pl -r ${OUTDIR}/outputDir/accepted_hits.sam -p ${OFFSET_FILE_PATH} -o ${OUTDIR}/outputDir/corrected.${SAMPLE}.mapping.sam
    perl ${RIBORF_CODE}/RibORF.2.0/readDist.pl -f ${OUTDIR}/outputDir/corrected.${SAMPLE}.mapping.sam -g ${HG19_GENEPRED} -o ${OUTDIR}/outputDir -d 1


    module load samtools
    samtools view -Sb ${OUTDIR}/outputDir/corrected.${SAMPLE}.mapping.sam > ${OUTDIR}/outputDir/corrected.${SAMPLE}.mapping.bam
    samtools sort ${OUTDIR}/outputDir/corrected.${SAMPLE}.mapping.bam -o ${OUTDIR}/outputDir/corrected.${SAMPLE}.mapping.sorted.bam
    samtools index ${OUTDIR}/outputDir/corrected.${SAMPLE}.mapping.sorted.bam

    samtools view -Sb ${OUTDIR}/outputDir/accepted_hits.sam > ${OUTDIR}/outputDir/accepted_hits.bam
    samtools sort ${OUTDIR}/outputDir/accepted_hits.bam -o ${OUTDIR}/outputDir/accepted_hits.sorted.bam
    samtools index ${OUTDIR}/outputDir/accepted_hits.sorted.bam

    perl ${RIBORF_CODE}/RibORF.2.0/ribORF.pl -f ${OUTDIR}/outputDir/corrected.${SAMPLE}.mapping.sam -c ${HG19_CANDIDATE_ORF} -o ${OUTDIR}/outputDir
    module purge

}


# offical tutorial: https://github.com/zhejilab/RibORF/tree/master/RibORF.2.0
# download software (riborf and gtfToGenePred) and necessary input (I use hg19.fa and hg19.ensGene.gtf from ucsc as I use them for nuORF extraction, if gencode, v19 and below for hg19)
# git clone https://github.com/zhejilab/RibORF/
# wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
# another material: https://docs.google.com/presentation/d/1cpkP9cXqOjiVIqoLpL8Hfw8PevqZUqNKX7DkDmGIGXI/edit?usp=sharing
# https://github.com/frankligy/NeoVerse/blob/main/archive/NeoVerse/internal/get_nuorf_fasta.py get search space from riborf

RIBORF_CODE=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/RibORF
GTF_TO_GENEPRED_PROGRAM=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/gtfToGenePred
HG19_SEQ=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/hg19.fa
HG19_GTF=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/hg19.ensGene.gtf
HG19_GENEPRED=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/hg19.ensGene.genePred.txt  # gtfToGenePred file.gtf file.genePred.txt
HG19_CANDIDATE_ORF=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/candidateORF.genepred.txt # perl ORFannotate.pl hg19.fa -t file.genePred.txt -o .
HG19_CANDIDATE_SEQ=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/candidateORF.fa  # same as above
FASTQ=/gpfs/data/yarmarkovichlab/Frank/pan_cancer/normal_ribo/raw/SRR15513149.fastq.gz
ROOT_DIR=/gpfs/data/yarmarkovichlab/Frank/pan_cancer/normal_ribo/result
ADAPTER_SEQ=AGATCGGAAG  # using fastqc sample.fastq.gz to infer, first 10nt
CUTADAPT=/gpfs/data/yarmarkovichlab/russell_data/cutadapt_env/bin/cutadapt
RRNA_FASTA=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/riborf_test/human.ribosomal.rna.fa
BOWTIE_HG19_INDEX_DIR=/gpfs/data/yarmarkovichlab/neuroblastoma/riboseq/ribo_real  # bowtie2-build file.fa hg19genome.index

SAMPLE=$(basename -s .fastq.gz ${FASTQ})
OUTDIR=${ROOT_DIR}/${SAMPLE}

if [ ! -d ${OUTDIR} ]; then
    mkdir ${OUTDIR}
fi

run_riborf_part1
# run_riborf_part2



