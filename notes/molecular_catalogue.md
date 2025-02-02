# Instructions on released molecular catalogue

We released the whole molecular catalogues for 21 cancers [here](https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/molecular_catalogue/), including all the modalties investigated in our study, including statitics on both sample and cohort level. We explain those files as below.

## Gene

1. `gene_tpm.txt`: tpm for each gene in each tumor sample.
2. `gene_lfc.txt`: cohort level gene expression for each gene and how tumor-specific they are compared to GTEx.

## Splicing

1. `splicing_all.txt`: record supporting read counts for each splicing junction in each tumor sample.
2. `splicing_rec.txt`: cohort level splicing statistics, how tumor specific they are, particularly we annotate the known splicing site and the associated genes, along with whether the splicing site fall into an annoated TE region (te_info column).

## Intron Retention

1. `intron_all.txt`: record intron expression (Stringtie reported exons coverage that read through an intron) in each tumor sample.
2. `intron_rec.txt`: cohort level and to what extent they are present in normal tissue.
3. `intron_peptide_all.txt`: all the read-through transcripts, their cds and peptide sequence.

## Mutation

1. `TCGA-cancer.mutect2_snv.tsv`: or similar suffix, documenting sample level somatic mutations
2. `mutation_rec.txt`: cohort level, how frequenty they are and whether they are driver gene mutation.

## Fusion

1. `fusion.txt`: sample level fusion presence.
2. `fusion_rec.txt`: cohort level fusion and whether they are detected in matched normal sample and other tumors, whether they are reported in literature or not, etc.

## Transposable Element

Note: Although the name contains erv, they actually contains all the TE elements. Sorry for the confusions.

1. `tumor_erv.h5ad`: This contains all the TE and the raw count and normalized cpm value in each tumor.

```python
import anndata as ad
tumor_erv = ad.read_h5ad('tumor_erv.h5ad')
tumor_erv.X.toarray() # will give you n_te * n_sample, raw count
tumor_erv.obs_names # will give you te names
tumor_erv.var_names # will give you tumor sample names
tumor_erv.layers['cpm'].to_array() # will give n_te * n_sample, cpm values
```

2. `ERV.txt`: cohort level, how tumor specific they are.
3. `good_erv.txt`: tumor specific TEs.

## HLA type

1. `hla_types`: sample level HLA-I information
