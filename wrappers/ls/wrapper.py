from snakemake.shell import shell
shell('ls {snakemake.input} > {snakemake.output}')
