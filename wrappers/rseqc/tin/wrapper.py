#!/usr/bin/env python

import os
from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')


# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
log = snakemake.log_fmt_shell()

shell(
    'tin.py '
    '-i {snakemake.input.bam} '
    '-r {snakemake.input.bed} '
    '{extra} '
    '{log}')

name = os.path.basename(snakemake.input.bam).rstrip('.bam')

shell(
        'mv {name}.tin.xls {snakemake.output.table} '
        '&& mv {name}.summary.txt {snakemake.output.summary}')
