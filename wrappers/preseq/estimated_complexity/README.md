# Wrapper for Preseq lc_extrap

Preseq method used to estimate the future yield of the genome library based on
observed data.

## Example

```python
rule preseq_lcextrap:
    input:
        bam='{sample}.sorted.bam'
    output:
        txt='{sample}/future_yield.txt'
    wrapper:
        "file://path/to/wrapper"
```

## Input

* `bam`: Coordinate-sorted BAM. If paired-end, include `-P` in params.extra to
  run in paired-end mode.

## Output

* `txt`: four column text file displaying the *total reads*, the
  *expected distinct* reads and its corresponding *lower/upper 95% Confidence
  Interval*

## Threads

Threads not supported.

## Params

* `extra`: passed verbatim; include `-P` to run in paired-end mode.
