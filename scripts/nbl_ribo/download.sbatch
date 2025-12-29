#!/bin/bash
#SBATCH --partition=a100_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=100Gb
#SBATCH --job-name="nbl_ribo"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:a100:1

./run_riborf.sh


# download ribo

# OHMX20210030_001.fastq.gz
# OHMX20210030_002.fastq.gz
# OHMX20210030_003.fasterq
# OHMX20210030_004.fasterq
# OHMX20210030_005.fastq.gz
# OHMX20210030_006.fastq.gz
# OHMX20210030_007.fastq.gz
# OHMX20210030_008.fasterq

# declare -a index_array=(OHMX20210030_001.fastq.gz OHMX20210030_002.fastq.gz OHMX20210030_003.fasterq OHMX20210030_004.fasterq \
#                         OHMX20210030_005.fastq.gz OHMX20210030_006.fastq.gz OHMX20210030_007.fastq.gz OHMX20210030_008.fasterq)
# for item in ${index_array[@]}
# do 
#     sshpass -p 'forPOSTDOC1!' sftp "guangyuan.li@nyulangone.org@transfer2.chop.edu:/NBL_Riboseq\ (farrela)/OHMX_Data/uploads/ribo/fastq/${item}" ./ribo
# done


# download rna

# OHMX20210030_001    
# OHMX20210030_002    
# OHMX20210030_003    
# OHMX20210030_004    
# OHMX20210030_005    
# OHMX20210030_006    
# OHMX20210030_007  
# OHMX20210030_008  

# declare -a index_array=(OHMX20210030_001 OHMX20210030_002 OHMX20210030_003 OHMX20210030_004 \
#                         OHMX20210030_005 OHMX20210030_006 OHMX20210030_007 OHMX20210030_008)
# for item in ${index_array[@]}
# do
#     sshpass -p 'forPOSTDOC1!' sftp "guangyuan.li@nyulangone.org@transfer2.chop.edu:/NBL_Riboseq\ (farrela)/OHMX_Data/uploads/rna/fastq/${item}/*" ./rna
# done