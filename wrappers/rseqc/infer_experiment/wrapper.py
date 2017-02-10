#!/usr/bin/env python

from tempfile import NamedTemporaryFile
from snakemake.shell import shell

# All wrappers must be able to handle an optional params.extra.
extra = snakemake.params.get('extra', '')


# This lets us handle whether to write to a log file or to write to stdout.
# See snakemake.script.log_fmt_shell for details.
log = snakemake.log_fmt_shell()

tmp = NamedTemporaryFile().name

shell(
    'infer_experiment.py '
    '-i {snakemake.input.bam} '
    '-r {snakemake.input.bed} '
    '{extra} '
    '> {tmp} '
    '{log}')

shell(
    'mv {tmp} {snakemake.output.txt}'
        )
