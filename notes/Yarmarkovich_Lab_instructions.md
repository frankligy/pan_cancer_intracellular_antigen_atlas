# Instructions of running ImmunoVerse for yarmarkovich lab members

## Before you start

This tutorial assumes you already have a bigpurple account, and have permission to access yarmarkovich lab drive, and have some basic understanding on how to run linux command and submit jobs through slurm scheduler. If you are not familiar with this, please carefully read through [yarmarkovich lab HPC tutorial](https://docs.google.com/document/d/14PZOJr-cr6z3JhHkEptWw5fPAztqVsmYXn3FvTKL4Zs/edit?usp=sharing).

This tutorial also assumes you have basic understanding of RNA-Seq data, immunopeptidome data, for instance, what is FASTQ file, what is read alignment step. Although following the steps should be straightford, I'd recommend you to use `chatGPT` for all any inquiry you may have as you go through this tutorial.

You should have some tumor RNA-Seq and some tumor immunopeptidome data at hand when you are reading this tutorial.

## Step 1: Organize your input data, codes, and outputs

I recommend you to have a folder structure like following:

```bash
/path/to/RNA_raw_data
    sample1_R1.fastq.gz
    sample1_R2.fastq.gz
    sample2_R1.fastq.gz
    sample2_R2.fastq.gz
/path/to/immunopeptidome_raw_data
    ./sample1
        sample1_rep1.raw
        sample1_rep2.raw
        sample1_rep3.raw
    ./sample2
        sample2_rep1.d
        sample2_rep2.d
        sample2_rep3.d
```

Besides the inputs, I also recommend you have the following folder, I will explain some of the files in a bit:

```bash
/path/to/codes
    template.json
    sub_all.sbatch
/path/to/result
    ./immunoverse_result
    samples.txt
```

The `template.json` and `sub_all.sbatch` are the keys for running the program, `json` file is an organized plain text file that specifity your input/output path, parameters for the program, this tutorial will primarily focus on how to modify this file. `sbatch` file contains the script to run the program and submit to the background.

I have created two template files, please copy them to your `codes` directory:

```bash
cp /gpfs/data/yarmarkovichlab/softwares/NeoVerse/template.json /path/to/codes
cp /gpfs/data/yarmarkovichlab/softwares/NeoVerse/sub_all.sbatch /path/to/codes
```

The `result` folder will be the place where the output will go into, within which there will be another subfolder called `immunoverse_result` which will contain more summarizied data, normally, all the results you care about should be in `immunoverse_result` folder. You don't have to manually create these two folders, but I want to illustrate in this way so you won't get too confused as you continue reading this tutorials. Another thing is to manually create a `samples.txt` file, it should be as following in this dummy case:

```bash
# samples.txt
sample1
sample2
```

## Step 2: Modify json file and run NeoVerse.sh

Simply put, this step involves running `NeoVerse.sh` script like this, these code block should be in your `sbatch` file, you just need to make them active and submit the job:

```bash
/gpfs/data/yarmarkovichlab/softwares/NeoVerse/NeoVerse.sh template.json
```

But you need to modify your `template.json` file to tell the program where your inputs are, where to generate your output, and how to run the program. Please locate to `shell_script` and `tunable_parameters` part, and change accordingly:

```json
    "intdir":"/path/to/RNA_raw_data",
    "rootdir":"/path/to/result",
    "sample":"sample1",
    "strand":"no",
    "mode":"normal_run",
    "modality":"rna"
```

One thing you may notice is here I set `sample` as `sample1`, that's because every single time you can only run one sample, so here you should submit one job for each sample, you should have 2 jobs running in parallel. 

For this step, please set `--cpus-per-task=3` and `--time=1-00:00:00` and `--mem=150Gb`.

After this step, You should have two folders generated under your `result` folder.

## Step 3: Modify json file and run NeoVerse.py 

Again, simply put, this step involves running `NeoVerse.py` script like below:

```bash
module load samtools
module load htslib/1.3
module load singularity/3.9.8
export LD_LIBRARY_PATH=/gpfs/data/yarmarkovichlab/Frank/immunopeptidome_project/engine/SNV/snv_env/lib:$LD_LIBRARY_PATH
/gpfs/data/yarmarkovichlab/softwares/NeoVerse/NeoVerse.py --config template.json --running_mode non_interactive
```

And again you need to modify your `json` file, please locate to `python_script` and `tunable_parameters`:

```json
    "min_pep_len":8,
    "max_pep_len":15,
    "outdir":"/path/to/result/immunoverse_result",
    "rootdir":"/path/to/result",
    "sample_file":"/path/to/result/samples.txt",
    "library":"no",
    "home_dir":"/path/to/your/homedir",
    "db_fasta_dir":"db_fasta",
    "image_format":"pdf",
    "kallisto_rootdir":"na",
    "kallisto_sample":"na"
```

After this step, you should have a folder called `db_fasta` in your `immunoverse_result`, which contains sample-specific search space, please use this to search your immunopeptidome data using Tesorai.

## Step 4: Generate final tabular output from Tesorai result

Placeholder

## Appendix I: How to interpret the source aberration name?

- Self gene

```bash
placeholder
```

- Splicing, TE chimeric transcript, self translated TE

```bash
# splicing junction (format 1): Junction coordinate (hs1 reference)|strand|Supporting RNA read count|phase|splicing site 1/2 gene symbols|nth discontinous peptide fragment from this junction
chr9:33946622-33951162|-|14|0|None,UBAP2|1

# splicing junction (format 2): Junction coordinate (hs1 reference)|strand|Supporting RNA read count|This junction is in a documented isoform ENST ID|splicing site 1/2 gene symbols
chr1:9774594-9775848|+|12|ENST00000377083|KIF1B,KIF1B

# TE chimeric transcript: Junction coordinate (hs1 reference)|strand|Supporting RNA read count|phase|TE element contributing to exonization (TE gene ID, TE transcript ID, TE family ID, TE class ID, server as splicing donor or acceptor, strand, coordinates)|splicing site 1/2 gene symbols|nth discontinous peptide
chr12:96898675-96908677|+|7|1|TE_info:L1ME2,L1ME2_dup1895,L1,LINE,donor,-,96897569,96898675|None,NEDD1|29

# self trasnlated TE: TE transcript ID|tumor read count in this sample|coordiante and strand|median tumor read count in tumor cohort|mean normal read count in gtex|fold change|# samples having this TE in tumor cohort|nth discontinous peptide fragment from this TE|phase|sense or anti-sense translation for TE
HERVH-int_dup2802|349.44|chr7:26026526-26029368:-|3.75|0.04|7781.16|4|71|3|sense
```

To explain `phase` in the context of splicing:

```bash
# phase 0
* * * * * *, * * * * * *
----- -----  ----- -----
# phase 1
* * * * * *, * * * * * *
        ------
# phase 2
* * * * * *, * * * * * *
          ------
```

To explain `phase` in self translated TE:

```bash
# phase 1, starting from first nt
* * * * * *
-----
# phase 2, starting from second nt
* * * * * *
  -----
# phase 3, starting from third nt
* * * * * *
    -----
```

- nuORF

```bash
# we adopted nuORF database convention, hg19 coordinate and the ENST transcript which carries this cryptic ORF
ENST00000534336.1_1_11:65268195-65268336:+|nuORF
```

- variant 

```bash
# ENSG ID|ENST ID|gene symbol|TPM expresion|mutation|type of mutation|likely RNA edit event or not|variant allele frequency or clonalities
ENSG00000206503|ENST00000376809|HLA-A|tpm:294.778174|L10V|missense_variant|possible_rna_edit_False|0.987297554779295
```







