__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()

if 'R1' in snakemake.input.keys() and 'R2' in snakemake.input.keys():
    shell(
        "cutadapt "
        "{extra} "
        "{snakemake.input.R2} "
        "{snakemake.input.R2} "
        "-o {snakemake.output.R1} "
        "-p {snakemake.output.R2} "
        "{log}"
    )
else:
    shell(
        "cutadapt "
        "{extra} "
        "{snakemake.input.fastq} "
        "-o {snakemake.output.fastq} "
        "{log}"
    )
