#!/usr/bin/env python

from tempfile import NamedTemporaryFile
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

shell(
    'geneBody_coverage.py '
    '-i {bam} '
    '-o tmp '
    '-r {snakemake.input.bed} '
    '{extra} '
    '{log}'
    )

shell(
    'mv tmp.geneBodyCoverage.r {snakemake.output.r} '
    '&& mv tmp.geneBodyCoverage.txt {snakemake.output.txt} '
    '&& mv tmp.geneBodyCoverage.curves.* {snakemake.output.img}'
    )
