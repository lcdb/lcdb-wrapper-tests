# Wrapper for `preseq c_curve`

Preseq method that generates a complexity plot of the genome library based on observed data.

## Examples

Single-end mode

```python
rule preseq_ccurve:
    input:
        bam='{sample}.sorted.bam'
    output:
        table='{sample}.complexity_output.txt'
    wrapper:
        "file://path/to/wrapper"
```

## Input
* `bam`: Coordinate-sorted BAM. If paired-end, include `-P` in params.extra to
  run in paired-end mode.

## Output
* `table`: Two-column text file displaying the *total reads* and
  a corresponding number of *distinct reads*

## Threads
Threads not supported.

## Params
* `extra`: passed verbatim; include `-P` to run in paired-end mode.
