# Overview

In case users want to run our CGC docker locally, we expose the path of entrypoint shell script, and how the inputs, arguments and paramters were passed to the shell script.

Noting currently we are not actively maintaining this approach, for experienced bioinformaticiane, they can dig into the shell script in docker container and make necessary modifications.

## Gene_pipeline

```bash
/usr/lib/gene.sh 3 . /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/kallisto_index no /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/ensembl_protein.fasta /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/uniprot_reviewed_curated_addition.fasta /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/nuorf.fasta /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_R1.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_R2.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN20-9844_R1.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN20-9844_R2.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_R1.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_R2.fastq.gz
```

## Alignment_pipeline

```bash
/usr/lib/alignment.sh 3 . /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/star_hg38_index /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/hg38.fa /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_R1.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_R2.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN20-9844_R1.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN20-9844_R2.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_R1.fastq.gz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_R2.fastq.gz
```

## Splicing_intron_pipeline

```bash
/usr/lib/splicing_intron_pipeline.sh 40 . /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/gene_model.txt /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/hg38.fa /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/gencode.v36.annotation.gtf no /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_secondAligned.sortedByCoord.out.bam /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN20-9844_secondAligned.sortedByCoord.out.bam /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_secondAligned.sortedByCoord.out.bam
```

## Pathogen_pipeline

```bash
/usr/lib/pathogen.sh 1 . /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/kraken2_db pair /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_secondAligned.sortedByCoord.out.bam /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN20-9844_secondAligned.sortedByCoord.out.bam /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_secondAligned.sortedByCoord.out.bam
```

## TE_local_pipeline

```bash
/usr/lib/te_pipeline.sh 20 . /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/hg38.knownGene.gtf /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/hg38_rmsk_TE.gtf.locInd no /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_secondAligned.sortedByCoord.out.bam /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN20-9844_secondAligned.sortedByCoord.out.bam /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_secondAligned.sortedByCoord.out.bam
```

## STAR_fusion (built by CGC)

```bash
tar -xvf /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/GRCh38.d1.vd1.gencode.v36.annotation.star-fusion-1.12.0-CTAT-index-archive.tar -C . && STAR-Fusion --left_fq /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_R1.fastq.gz --right_fq /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_R2.fastq.gz --genome_lib_dir `pwd`/ctat_genome_lib_build_dir --output_dir fusion_output  --examine_coding_effect  --CPU 32   --tmpdir /tmp --verbose_level 2 && tar -vcf HN21.fusion_output.tar fusion_output && mv fusion_output/star-fusion.fusion_predictions.abridged.tsv fusion_output/HN21.star-fusion.fusion_predictions.abridged.tsv && mv fusion_output/star-fusion.fusion_predictions.tsv fusion_output/HN21.star-fusion.fusion_predictions.tsv && if [ -e fusion_output/FusionInspector-undefined/finspector.FusionInspector.fusions.abridged.tsv ]; then mv fusion_output/FusionInspector-undefined/finspector.FusionInspector.fusions.abridged.tsv fusion_output/FusionInspector-undefined/HN21.finspector.FusionInspector.fusions.abridged.tsv; fi && if [ -e fusion_output/FusionInspector-undefined/finspector.fusion_inspector_web.html ]; then mv fusion_output/FusionInspector-undefined/finspector.fusion_inspector_web.html fusion_output/FusionInspector-undefined/HN21.finspector.fusion_inspector_web.html; fi
```

## Variant_pipeline (Step 1: Get VCF from BAM)

```bash
/usr/lib/variant_pipeline.sh 3 . /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/GRCh38.d1.vd1.fa pair /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_secondAligned.sortedByCoord.out.bam /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN20-9844_secondAligned.sortedByCoord.out.bam /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_secondAligned.sortedByCoord.out.bam
```

## Variant_pipline (Step 2: Run VEP to get coding effect, built by CGC)

```bash
tar xfz /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/homo_sapiens_vep_112_GRCh38.tar.gz -C /opt/ensembl-vep/ && perl -I /root/.vep/Plugins/ /opt/ensembl-vep/vep --buffer_size 20000 --fork 8 --input_file /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN21-10181_variants.vcf --offline --dir /opt/ensembl-vep --output_file HN21-10181_variants.vep.vcf --vcf ; sed -i 's=http://www.google.com/jsapi=https://www.google.com/jsapi=g' *summary.html
```

## HLA typing using Optitype (built by CGC)

```bash
# skip first step which convert fastq.gz to fastq
python /opt/OptiType-1.3.5/OptiTypePipeline2.py --input /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_R1.fastq /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/HN19-9674_R2.fastq --rna --outdir ./temp/ --config ./config.ini > HN19.command_log.txt && mv ./temp/*/*_result_type.tsv HN19.result_type.tsv && mv ./temp/*/*_coverage_plot.pdf HN19.coverage_plot.pdf && mv ./temp/*/*_result.tsv HN19.result.tsv && mv ./temp/*/*_result_id.tsv HN19.result_id.tsv && mv ./config.ini HN19_config.ini
```

## Summarization Step

```bash
/usr/lib/summary.sh /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/result /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/ImmunoVerse_data .
```

## Immunopeptidome 1: MaxQuant

```bash
/usr/lib/run.sh /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/test_immuno /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/test_fasta orbitrap 1 1 . hg19
```

## Immunopeptidome 2: MSconvert

```bash
wine msconvert /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/test_immuno/20240110_E_OdinLC_IC_PDX_HD_19.raw --outdir . --zlib --filter "peakPicking vendor msLevel=1-"
```

## Immunopeptidome 3: Rescore

```bash
/usr/lib/rescore.sh /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/test_mzml /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/hg19_maxquant_combined_txt . 0.05 orbitrap test
```

## Immunopeptidome 4: netMHCpan binding prediction

```bash
/usr/lib/binding.sh /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/test_msmsScans_new.txt /sbgenomics/Projects/fc4dcac6-0e95-4c65-8558-0b8d559353f1/test_hla_type.txt . 1 test_final2
```