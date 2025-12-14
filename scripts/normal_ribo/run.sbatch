#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=100Gb
#SBATCH --job-name="ribo"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err


# module load sratoolkit/2.10.7

# # GSE182377 n=79
# # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5527673 AGATCGGAAG
# grep "ribosome-protected mRNA fragments" SraRunTable_GSE182377.csv | cut -d "," -f 1 | xargs -L 1 -P 10 -I {} sh -c "fasterq-dump -e 20 {}"
# for file in *.fastq; do echo $file; done | xargs -L 1 -P 10 -I {} sh -c "gzip {}"

# # GSE59820 n=4
# # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1446874 TGGAATTCTC
# declare -a samples=(
# SRR1528686
# SRR1528687
# SRR1528688
# SRR1528689
# )

# for s in ${samples[@]}; do echo $s; done | xargs -L 1 -P 10 -I {} sh -c "fasterq-dump -e 20 {}"
# for file in *.fastq; do echo $file; done | xargs -L 1 -P 10 -I {} sh -c "gzip {}"

./run_riborf.sh
