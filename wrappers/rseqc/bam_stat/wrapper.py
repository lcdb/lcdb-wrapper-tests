#!/usr/bin/env python

from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')

# bam_stat prints the desired output to STDERR so don't log

shell(
    'bam_stat.py '
    '-i {snakemake.input.bam} '
    '{extra} '
    '2> {snakemake.output.txt} ')
