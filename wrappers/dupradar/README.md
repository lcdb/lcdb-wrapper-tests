# Wrapper for dupRadar

dupRadar provides an easy way to distinguish between artifactual vs natural
duplicate reads in RNA-Seq data. Prior to dupRadar only global duplication rates
were used and they don't take into account the effect of gene expression levels. 
dupRadar relates *duplication rates* and *length normalized read counts* of every
gene to model the dependency of both variables. 

[Link to homepage](https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html)

[Link to manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html)

## Example

Single-end, not stranded:

```python
rule dupRadar:
    input:
       bam='{sample_dir}/{sample}/{sample}.cutadapt.hisat2.unique.sort.dedup.bam',
       annotation='annotations/dm6.gtf',
    output:
        png='{sample_dir}/{sample}/{sample}_dupRadar_drescatter.png'
    wrapper:
        wrapper_for('dupRadar')
```

Paired-end, stranded:

```python
rule dupRadar:
    input:
       bam='{sample_dir}/{sample}/{sample}.cutadapt.hisat2.unique.sort.dedup.bam',
       annotation='annotations/dm6.gtf',
    output:
        png='{sample_dir}/{sample}/{sample}_dupRadar_drescatter.png'
    params:
        paired=True,
        stranded=True
    wrapper:
        wrapper_for('dupRadar')
```

## Input
* `bam`: BAM file with mapped reads has to be duplicate marked using either
  Picard or BamUtil

* `annotation`: GTF file contaning features to count the reads falling on the
  features.

## Output
* `png`: an expression density scatter plot in .png format is generated

## Threads
Threads are passed to dupRadar.

## Params
* `paired`: True | False. Default False.
* `stranded`: True | False | "reverse". Default False.
