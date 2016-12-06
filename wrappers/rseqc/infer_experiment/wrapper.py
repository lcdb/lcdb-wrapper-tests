#!/usr/bin/env python

from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')


# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
log = snakemake.log_fmt_shell()


shell(
    'infer_experiment '
    '-i {snakemake.input.bam} '
    '-r {snakemake.input.bed} '
    '{snakemake.params.extra} '
    '> {snakemake.output[0]} '
    '{log}')
