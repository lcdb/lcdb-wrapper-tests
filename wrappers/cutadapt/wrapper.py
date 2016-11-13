__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

from snakemake.shell import shell
try:
    extra = snakemake.params.extra
except AttributeError:
    extra = ""
log = log_fmt_shell()
shell(
    "cutadapt "
    "{extra} "
    "{snakemake.input.fastq} "
    "-o {snakemake.output.fastq} "
    "{log}"
)
