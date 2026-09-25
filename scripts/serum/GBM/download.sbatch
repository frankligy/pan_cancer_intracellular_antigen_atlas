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



# PXD008127_plasma_and_tumor just plasma
declare -a samples=(
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq35945_QE.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39926_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq35947_QE.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39927_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39929_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43508_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39931_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43509_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43482_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43481_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43484_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43483_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43486_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43485_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43490_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43489_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39924_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39925_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43491_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43493_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43492_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43495_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43494_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43464_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43463_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43496_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47972_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47974_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47977_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47940_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47937_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47934_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47936_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47935_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47939_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47938_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47933_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47948_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47952_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47955_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47951_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47954_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47950_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47953_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47949_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47957_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47956_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47964_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47960_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47963_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47959_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47962_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47958_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47961_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47966_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47965_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47947_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47943_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47946_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47945_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47944_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47942_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47976_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43475_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47971_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47970_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43504_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43505_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47969_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47968_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq47967_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43477_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43506_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq35946_QE.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq35949_QE.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39928_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq35948_QE.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39930_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43507_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39933_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39934_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq39935_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43510_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43511_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43512_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43480_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43479_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43458_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43457_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43488_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43487_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43461_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43460_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43467_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43466_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43498_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43497_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43499_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43501_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43500_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43470_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43469_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43473_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43472_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43503_QE2.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/06/PXD008127/Seq43502_QE2.raw
)

for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P ../immuno/dummy {}"


