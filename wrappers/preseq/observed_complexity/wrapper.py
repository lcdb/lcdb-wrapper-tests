#!/usr/bin/env python

from snakemake.shell import shell
log = snakemake.log_fmt_shell()
extra = snakemake.params.get('extra', '')

shell("preseq c_curve {extra} -B {snakemake.input.bam} -o {snakemake.output.table}")
