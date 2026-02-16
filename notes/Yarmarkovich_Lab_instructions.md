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

**Note: This section is more for myself, I don't expect you to use the following code, instead, please contact Ritchlynn, Xinya or Aman on how to do the post-hoc analysis, mainly how to dissect differnet types of antigens, and how to conduct binding predictions.**

Again, the actual code you need to run is as easy as:

```bash
/gpfs/data/yarmarkovichlab/softwares/NeoVerse/post_maxquant.py --config template.json
```

Hard part is to modify the `json` file, please locate `maxquant` and `tunable_parameter` section, it's going to intimating at first glance, but most of the parameters can be left untouched (the one I left as na). And for now, I assume we will be using the tesorai `quantified_psm` file, you shall remove the file extension (so tesorai will report file.zip or file.raw, please remove the file extension using excel or writing code) under the column `filename` (no underscore), and put this processed file as `/path/to/immunopeptidome_raw_data/sample_x/combined/txt/other_alg.txt`, my program will particularly look for this path to pick up the result file.

```json
    # technology should be keyed on rna sample name
    "outdir":"/path/to/result/immunoverse_result",
    "db_fasta":"db_fasta",
    "immunopeptidome_dir":"/path/to/immunopeptidome_raw_data",
    "technology":{
        "rna_sample1":"orbitrap",
        "rna_sample2":"bruker"
    },
    "rna2immuno":{
        "rna_sample1":"immuno_sample1",
        "rna_sample2":"immuno_sample2"
    },
    "peptide_fdr":"na",
    "hla_class":"na",
    "dia_dir":"na",
    "tmt_info":"na",

    "current_maxquant_run_mqpar":"na",


    "rederive_mzml_dir":"na",
    "rederive_vanilla_fdr":"na",
    "rederive_mode":"na",
    "rederive_train_cutoff":"na",
    "rederive_rescore_fdr":"na",
    "rederive_samples":"na",

    "mode":"other_alg",
    "hla_mapping":{
        "rna_type_sample1":"immuno_sample1",
        "rna_type_sample2":"immuno_sample2"
    },
    "overwrite_hla_dic":null,
    "overwrite_hla_nested":false,
    "overwrite_additional_hla":false,
    "additional_hla":[],
    "cores":20,
    "inquiry_mode":"i",
    "added_genes":[],
    "use_genes_lfc":true, # if not rna, use false
    "use_bayesTS":true, # if not rna, use false
    "intensity_as_dict":false,
    "other_alg_mapping":{
        "filename":"Raw file",
        "clean_sequence":"Sequence",
        "scan_id":"Scan number",
        "protein_ids":"Proteins",
        "intensity":"Precursor intensity",
        "score":"Score",
        "qval":"PEP",
        "precursor_charge":"Charge",
        "precursor_mz":"m/z",
        "retention_time":"Retention time"
    },
    "other_alg_impute":{
        "plot_logic":"mzml",
        "Mass analyzer":"Bruker_TIMS_TOF"
    },
    "protein_delimiter":";"

```

In certain cases, the correspondance between your RNA and immunopeptidome is not exactly one-to-one or even more complexed. You may have noticed `overwrite_hla_dic` and `overwrite_hla_nested` parameters to allow you to specify HLA types for each immunopeptidome samples instead of reading from your RNA-Seq based pipeline. You should set `hla_mapping` to `null` as you don't need this mapping anymore.

* Scenario 1:

```json
# [] means no HLA annotaiton for that immunopeptidome sample
"overwrite_hla_dic":{
    "immuno_sample1":["A*32:01","B*40:01","C*03:04"],
    "immuno_sample2":[]
}
```

* Scenario 2:

```json
# Sometimes, we group multiple samples in one study
"overwrite_hla_dic":{
    "immuno_sample1@actual_sample1.raw":["A*32:01","B*40:01","C*03:04"],
    "immuno_sample1@actual_sample2.raw":["A*31:01","A*32:01","B*14:01","B*27:05","C*02:02","C*08:02"],
    "immuno_sample2@actual_sample1":[]
},
"overwrite_hla_nested":true,
```

## Step 5: Visualization

**Notes: This part is written for me, if you'd like to visualize the PSM, consult tesorai team, if you'd like to visualize the differential plots, consult Aman**

Again, the actual code you need to run is as easy as:

```bash
/gpfs/data/yarmarkovichlab/softwares/NeoVerse/launch_portal.py --config template.json --running_mode generate_figures
```

Make sure you requested mzml from tesorai, and put them into `/path/to/tesorai_mzml`, folder structure should match immunopeptidome raw data folder, file extension should be `mzML`, then modify the json file as follow.

```json
    "raw_dir":"/path/to/immunopeptidome_raw_data",
    "mzml_dir":"/path/to/tesorai_mzml",
    "technology_setup":"bruker",
    "cores":20,

    "cancer":"cancer_alias",
    "antigen_dir":"/path/to/result/immunoverse_result/antigen/other_alg",
    "assets_dir":"/path/to/result/immunoverse_result/assets",

    "template_json":"/path/to/codes/template.json",
    "draw_diff":false,
    "draw_psm":true

```

And remember, you have to modify the `python_interact/nuorf` part, you can leave others unchanged, but the `raw2bio`, see below.

```json
    "nuorf":{
        "obj":"SLFEGIYTI",

        "immuno_dir":"/gpfs/data/yarmarkovichlab/JH_AML/immuno",
        "col":"percentile",

        "raw2bio":{
            "each_raw_file.raw":"immuno_sample1",
            "each_raw_file.d":"immuno_sample2"
        },

        "final_path":"/gpfs/data/yarmarkovichlab/JH_AML/antigen/other_alg/final_enhanced.txt"

    }
```

Lastly, you may overwrite `all_peps` and `png` in the `launch_portal.py` code temporarily, but also revert changes once you are done.


## Appendix I: How to interpret the source aberration name?

- Self gene

```bash
# ENSG|ENST|gene_symbol|TPM in this sample
ENSG00000138182|ENST00000371728|KIF20B|tpm:29.147336
```

- Splicing, TE chimeric transcript, self translated TE, circular RNA

```bash
# splicing junction (format 1): Junction coordinate (hs1 reference)|strand|Supporting RNA read count|phase|splicing site 1/2 gene symbols|nth discontinous peptide fragment from this junction
chr9:33946622-33951162|-|14|0|None,UBAP2|1

# splicing junction (format 2): Junction coordinate (hs1 reference)|strand|Supporting RNA read count|This junction is in a documented isoform ENST ID|splicing site 1/2 gene symbols
chr1:9774594-9775848|+|12|ENST00000377083|KIF1B,KIF1B

# TE chimeric transcript: Junction coordinate (hs1 reference)|strand|Supporting RNA read count|phase|TE element contributing to exonization (TE gene ID, TE transcript ID, TE family ID, TE class ID, server as splicing donor or acceptor, strand, coordinates)|splicing site 1/2 gene symbols|nth discontinous peptide
chr12:96898675-96908677|+|7|1|TE_info:L1ME2,L1ME2_dup1895,L1,LINE,donor,-,96897569,96898675|None,NEDD1|29

# self trasnlated TE: TE transcript ID|tumor read count in this sample|coordiante (hg38) and strand|median tumor read count in tumor cohort|mean normal read count in gtex|fold change|# samples having this TE in tumor cohort|nth discontinous peptide fragment from this TE|phase|sense or anti-sense translation for TE
HERVH-int_dup2802|349.44|chr7:26026526-26029368:-|3.75|0.04|7781.16|4|71|3|sense

# circular rna: Junction coordinate (hg38)|strand|phase|rna_count|nth discontinous peptide|circRNA
chr11:65499306-65500372|+|0|1|seq:2|circRNA
```

To explain `phase` in the context of splicing and circular RNA, `,` means junction:

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







