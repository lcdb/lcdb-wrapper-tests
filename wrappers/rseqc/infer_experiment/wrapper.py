#!/usr/bin/env python

from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')


# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
log = snakemake.log_fmt_shell()


# This demo shows how to handle paired-end and single-end input data as two
# different cases, depending on whether the rule's input included an "R2" key
# or not.

shell('infer_experiment -r {input.bed} -i {input.bam} {params.extra} > {output[0]}')
