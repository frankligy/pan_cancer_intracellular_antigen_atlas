# Instructions on released molecular catalogue

We also released the normal controls we used in the study [here](https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/normal/). Due to the massive space
for all GTEx sample (running on all sample for one modality can take months and also the storage charge is not always amenable for an academic lab), a normal practice is to randomly subsample equal amount of samples from each histology/tissue as a representative cohort. 

## Gene

* File: `gtex_gene_all.h5ad`
* Sample Size: GTEx (n=17,382)

```python
import anndata as ad
adata = ad.read_h5ad('gtex_gene_all.h5ad')  # n_gene × n_sample
adata.obs_names = [item.split('.')[0] for item in adata.obs_names]
adata.obs_names_make_unique()
```

## Splicing and Intron (Based on [SNAF](https://github.com/frankligy/SNAF) Database)

* File: `GTEx_junction_counts.h5ad`, `tcga_matched_control_junction_count.h5ad`
* Sample Size: GTEx (n=2,638), TCGA matched controls (n=701), most of GTEx tissues have around 25 samples
* The h5ad file use AltAnalyze internal identifiers, we also provided the correspoinding hg38 coordinates in bed format ([GTEx](https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/normal/gtex_tmp_prelift.bed), [TCGA](https://genome.med.nyu.edu/public/yarmarkovichlab/ImmunoVerse/normal/tcga_tmp_prelift.bed)), REMEMBER, bed format the start position is 0-based.

```python
# both are n_splicing * n_sample
adata_gtex = ad.read_h5ad('GTEx_junction_counts.h5ad')
adata_tcga = ad.read_h5ad('tcga_matched_control_junction_count.h5ad')
```

## Intron (Based on This paper, using StringTie)

* File: `normal_intron.txt`
* Sample Size: GTEx (n=155), each GTEx histology has 5 samples


## Pathogen

* File: `normal_pathogen.txt`
* Sample Size: GTEx (n=155), each GTEx histology has 5 samples

## Transposable Element

* File: `normal_erv.txt`
* Sample Size: GTEx (n=155), each GTEx histology has 5 samples

