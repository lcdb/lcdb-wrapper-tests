#!/usr/bin/env python

from tempfile import mkdtemp
from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')


# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
log = snakemake.log_fmt_shell()

# Make bam string
if isinstance(snakemake.input.bam, str):
    bam = snakemake.input.bam
elif hasattr(snakemake.input.bam, __iter__):
    bam = ','.join(snakemake.input.bam)
else:
    raise ValueError("BAM must be a single file or a list of files.")

# Make temp prefix
temp = mkdtemp


shell(
    'geneBody_coverage.py '
    '-i {bam} '
    '-o {temp} '
    '-r {snakemake.input.bed} '
    '{snakemake.params.extra} '
    '{log}')

shell(
    ''
    '-o {snakemake.output[0]} '
    '')
