#!/bin/bash
#SBATCH --partition=a100_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=150Gb
#SBATCH --job-name="PRIME"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:a100:1


# MixMHCpred

# conda activate -p ./MixMHCpred_env python=3.10
# ./MixMHCpred_env/bin/python3.10 -m pip install --no-cache-dir numpy pandas scipy argparse logomaker matplotlib
# conda install bioconda::mafft
# donwload the https://github.com/GfellerLab/MixMHCpred/releases/tag/v3.0 and unzip it, executable MixMHCpred

# PRIME

# download https://github.com/GfellerLab/PRIME/releases/tag/V2.1 and unzip it, executable PRIME
# cd ./PRIME-2.1/lib
# g++ -O3 PRIME.cc -o PRIME.x


module load anaconda3
conda activate /gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/webtool/PRIME/MixMHCpred_env
/gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/webtool/PRIME/PRIME-2.1/PRIME -i ./PRIME_input.txt -o ./PRIME_output.txt \
-a B4002,B4803,C0404,B5501,A6801,C0704,B3801,A0203,B3701,B2713,A6601,B2706,A3401,B1501,B3901,B4801,B5101,A2902,C0112,B5401,B3505,C0401,A0217,A2901,C0501,B3906,C0202,C0102,B5701,C1601,B4701,B5201,C0718,A0103,C0602,A0265,C1202,A7401,A3303,B1801,B4501,A0301,B5302,B5001,B3802,C07133,A2403,B2705,B4102,C0303,A02264,C1505,B5301,B1516,B1517,C1513,C0304,B4901,B1502,B5102,B4402,A3101,A0207,C1203,B0710,B5601,B3508,C1502,B0702,B3905,B5703,B3501,C0332,B1521,B4032,B4006,A3201,C0802,B3503,C1602,C0801,C0302,B0705,C0403,C1801,B1525,A2410,A0202,A0101,B1301,B7801,C1402,B3530,B2701,B4403,B5801,B4001,B2702,A3402,B4429,B5604,B1511,A0205,A0209,A0201,A2402,C0103,B35118,B1512,C1506,B4405,A1101,A3002,B0801,C0706,A3301,A2301,C1701,A3001,B4101,C0701,A0206,A2501,C0702,B3502,A2601,B1802,B1503,B1402,B4601,B1401,B5502,A6802,B1302 \
-mix /gpfs/data/yarmarkovichlab/Frank/pan_cancer/codes/webtool/PRIME/MixMHCpred-3.0/MixMHCpred