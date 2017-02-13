#!/usr/bin/env python

import os
from tempfile import NamedTemporaryFile
from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')

# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
log = snakemake.log_fmt_shell()

# tin uses the name of the BAM to create outputs. In order to write outputs to
# tmp I need to copy the BAM over to tmp.
bam = NamedTemporaryFile(suffix='.bam').name
bed = NamedTemporaryFile(suffix='.bed').name

shell(
    'cp {snakemake.input.bam} {bam} '
    '&& cp {snakemake.input.bam}.bai {bam}.bai '
    '&& cp {snakemake.input.bed} {bed}')

shell(
    'cd $TMPDIR '
    '&& tin.py '
    '-i {bam} '
    '-r {bed} '
    '{extra} '
    '{log}')

name = bam.rstrip('.bam')

shell(
        'mv {name}.tin.xls {snakemake.output.table} '
        '&& mv {name}.summary.txt {snakemake.output.summary}')
