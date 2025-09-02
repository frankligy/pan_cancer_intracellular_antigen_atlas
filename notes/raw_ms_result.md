# Instructions on released raw MS results

We released the [raw MS results (directly from searching algorithms)](https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/raw_MS_result/) to the community to enable reanalysis. Since diverse analytical procedures are often employed, including but not limited to FDR stringency, FDR calculation methods, etc, and their performances are often situation-dependent, we believe releasing the raw results would be the most proper and trasparent way to advance the field.

Specially, since we employed three procedures, namely `maxquant`, `MS2rescore` and `tesorai`, we explained the directory structure as below.

## Tesorai
Each cancer has one single tsv file. Most of the columns are self-explainable, `qval` is the FDR corrected q-value used to control FDR stringency.

## MaxQuant
We ran MaxQuant by batches for parallelization which are more scalable. All the batches information are available in our `Table S3`. Each batch has one dedicated subfolder, two main tables, namely `msms.txt` and `msmsScans.txt` are included for each run. We refer readers to the [offical website](https://cox-labs.github.io/coxdocs/output_tables.html) on each column in the table.

## MS2Rescore
We ran MS2Rescore on top of each MaxQuant batch, `ms2rescore.mokapot.psms.txt` file was included in each rescore run. `mokapot q-value` was used for FDR estimation as confirmed by the author (https://github.com/CompOmics/ms2rescore/discussions/196).

