# rseqc Infer Experiment

[website](http://rseqc.sourceforge.net/#infer-experiment-py)

## Examples

Minimal usage:

```python
rule infer_experiment:
    input: bam = 'a.bam',
           bed = 'b.bed'
    output: 'b.txt'
    log: 'a.infer_experiment.log'
    wrapper:
        'file://path/to/wrapper'
```

## Input

`bam`:
    Input alignment file in BAM format.

`bed`:
    Reference gene model in bed format.

## Output

Output files are simply copies of input.

### Single-end mode:

Expects a single unnamed output file

## Threads
Does not use threads

## Params
`extra`:
    Command line parameters that will be passed directly to rseqc infer_exerpiment.py

    Possible command include::
    
        -s SAMPLE_SIZE, --sample-size=SAMPLE_SIZE
            Number of reads sampled from BAM [200000].
        -q MAP_QUAL, --mapq=MAP_QUAL
            Minimum mapping quality (phred scaled) for an alignment to be considered [30].

