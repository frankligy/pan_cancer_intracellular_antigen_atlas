#!/bin/bash
#SBATCH --partition=a100_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=150Gb
#SBATCH --job-name="neoverse-pan"
#SBATCH --output=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.out
#SBATCH --error=/gpfs/data/yarmarkovichlab/Frank/job_dump/%j_%x.err
#SBATCH --gres=gpu:a100:1


# bladder
declare -a samples=(
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_1_400-650mz_msms32_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_2_400-650mz_msms33_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_3_400-650mz_msms34_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/161230_AM_AUT01-DN13_Bladder_W6-32_10_DDA_1_400-650mz_msms29_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/161230_AM_AUT01-DN13_Bladder_W6-32_10_DDA_2_400-650mz_msms30_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/161230_AM_AUT01-DN13_Bladder_W6-32_10_DDA_3_400-650mz_msms31_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180117_AM_AUT01-DN17_Bladder_W6-32_10_DDA_1_400-650mz_msms20_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180117_AM_AUT01-DN17_Bladder_W6-32_10_DDA_2_400-650mz_msms21_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180117_AM_AUT01-DN17_Bladder_W6-32_10_DDA_3_400-650mz_msms22_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/171109_AM_AUT01-DN11_Bladder_W6-32_10_DDA_1_400-650mz_msms1_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/171109_AM_AUT01-DN11_Bladder_W6-32_10_DDA_2_400-650mz_msms2_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/171109_AM_AUT01-DN11_Bladder_W6-32_10_DDA_3_400-650mz_msms3_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180317_AM_AUT01-DN09_Bladder_W6-32_10_DDA_1_400-650mz_msms17_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180317_AM_AUT01-DN09_Bladder_W6-32_10_DDA_2_400-650mz_msms18_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/180317_AM_AUT01-DN09_Bladder_W6-32_10_DDA_3_400-650mz_msms19_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170810_AM_AUT01-DN15_Bladder_W6-32_10_DDA_1_400-650mz_msms1_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170810_AM_AUT01-DN15_Bladder_W6-32_10_DDA_2_400-650mz_msms2_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170810_AM_AUT01-DN15_Bladder_W6-32_10_DDA_3_400-650mz_msms3_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170527_AM_AUT01-DN14_Bladder_W6-32_10_DDA_1_400-650mz_msms30_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170527_AM_AUT01-DN14_Bladder_W6-32_10_DDA_2_400-650mz_msms31_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170527_AM_AUT01-DN14_Bladder_W6-32_10_DDA_3_400-650mz_msms32_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170702_AM_AUT01-DN12_Bladder_W6-32_10_DDA_1_400-650mz_msms1_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170702_AM_AUT01-DN12_Bladder_W6-32_10_DDA_2_400-650mz_msms2_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/170702_AM_AUT01-DN12_Bladder_W6-32_10_DDA_3_400-650mz_msms3_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160927_AM_AUT01-DN06_Bladder_W6-32_10_DDA_1_400-650mz_msms13_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160927_AM_AUT01-DN06_Bladder_W6-32_10_DDA_2_400-650mz_msms14_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/160927_AM_AUT01-DN06_Bladder_W6-32_10_DDA_3_400-650mz_msms15_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/161218_AM_AUT01-DN08_Bladder_W6-32_10_DDA_1_400-650mz_msms15_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/161218_AM_AUT01-DN08_Bladder_W6-32_10_DDA_2_400-650mz_msms16_standard.raw
ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD019643/161218_AM_AUT01-DN08_Bladder_W6-32_10_DDA_3_400-650mz_msms17_standard.raw
)

for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/bladder {}"

# hepatocytes
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2064_L_MdB072_7mei19_donor_MdB026_OT.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2064_L_MdB072_7mei19_donor_MdB026_OT_2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2047_L_MdB_6jan2020_MdB024.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2047_L_MdB_6jan2020_MdB024_20200106173725.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2062_L_MdB070_3mei19_donor_MdB028_OT.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2062_L_MdB070_3mei19_donor_MdB028_2_OT.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2047_L_MdB_6jan2020_MdB031.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2047_L_MdB_6jan2020_MdB031_20200107003149.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2064_L_MdB072_7mei19_donor_MdB051_OT.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2064_L_MdB072_7mei19_donor_MdB051_OT_2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2062_L_MdB070_3mei19_donor_MdB052_OT.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/09/PXD023143/2062_L_MdB070_3mei19_donor_MdB052_2_OT.raw
)

for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/hepatocytes {}"

# beta_cell
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R1_Basal_even.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R1_Basal_odd.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R1_INFa_even.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R1_INFa_odd.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R2_Basal_even.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R2_Basal_odd.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R2_INFa_even.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R2_INFa_odd.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R3_Basal_even.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R3_Basal_odd.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R3_INFa_even.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R3_INFa_odd.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R4_Basal_even.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R4_Basal_odd.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R4_INFa_even.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD045265/FL1100_IPP132_R4_INFa_odd.raw
)

for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/beta_cell {}"

# iPSC
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/CB_iPSC_Ctl_251119_1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/CB_iPSC_Ctl_251119_2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/CB_iPSC_IFN_251119_1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/CB_iPSC_IFN_251119_2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/Fibro_iPSC_Ctl_251119_1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/Fibro_iPSC_Ctl_251119_2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/Fibro_iPSC_IFN_251119_1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/Fibro_iPSC_IFN_251119_2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/iPSC22_IFN_280619_1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/iPSC22_IFN_280619_2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/iPSC_375M_290119_1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/08/PXD026702/iPSC_375M_290119_2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/iPSC {}"


# Treg
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/03/PXD007860/FL0012779.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/03/PXD007860/FL0012783.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/03/PXD007860/FL0012803.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/03/PXD007860/FL0012807.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/03/PXD007860/FL0013327.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/03/PXD007860/FL0013331.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/03/PXD007860/FL0013353.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2018/03/PXD007860/FL0013357.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/Treg {}"

# CD14
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_CD14_3_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD14_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD14_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD14_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD14_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD14_3_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/CD14 {}"

# ImmDC
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_ImDC_3_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_ImmDc_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_ImmDc_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_ImmDc_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_ImmDc_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_ImmDC_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_ImmDC_1_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/ImmDC {}"

# MatureDC
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20171216_QEh1_LC1_HLAIp_SA_FaMa_Leuka2_MDC_3_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_MaDc_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_MaDc_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_MaDc_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_MaDc_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_MaDC_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_MaDC_1_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/MatureDC {}"

# CD4
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_02_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_02_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_03_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD4_03_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_02_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_02_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_03_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD4_03_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_3_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/CD4 {}"

# CD8
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_02_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_02_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD8_03_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD8_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD8_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD8_02_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD8_02_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_3_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/CD8 {}"

# CD19
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD19_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180424_QEh1_LC1_FaMa_HLAIp_D2_CD19_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD19_01_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD19_01_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD19_02_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180428_QEh1_LC1_FaMa_HLAIp_D3_CD19_02_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_3_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/CD19 {}"

# CD4_Act
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD4_Act_3_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/CD4_Act {}"

# CD8_Act
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_1_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_2_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_2_R2.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_3_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD8_Act_3_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/CD8_Act {}"

# CD19_Act
declare -a samples=(
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_Act_1_R1.raw
https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD020011/20180814_QEh1_LC1_SA_FaMa_HLAIp_CD19_Act_1_R2.raw
)
for item in ${samples[@]}; do echo ${item}; done | xargs -n 1 -L 1 -P 10 -I {} bash -c "wget -P /gpfs/data/yarmarkovichlab/Frank/pan_cancer/safety_screen/immuno/CD19_Act {}"