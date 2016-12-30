import tempfile
from snakemake.shell import shell
from lcdblib.snakemake import helpers

extra = snakemake.params.get('extra', '')
log = snakemake.log_fmt_shell()

stranded = snakemake.params.get('stranded', False)
try:
    stranded_int = {False: 0, True: 1, 'reverse': 2}[stranded]
except KeyError:
    raise ValueError('"stranded" must be True|False|"reverse"')

paired = snakemake.params.get('paired', False)
try:
    paired_bool= {True: 'TRUE', False: 'FALSE'}[paired]
except KeyError:
    raise ValueError('"paired" must be True or False')


script = """
library(dupRadar)
bam <- "{snakemake.input.bam}"
gtf <- "{snakemake.input.gtf}"
dm <- analyzeDuprates(bam, gtf, {stranded_int}, {paired_bool}, {snakemake.threads})
png(file="{snakemake.output.png}", width=1000, height=1000)
duprateExpDensPlot(dm, main=basename(bam))
dev.off()
""".format(**locals())

tmp = tempfile.NamedTemporaryFile(delete=False).name
helpers.rscript(script, tmp, log=log)
