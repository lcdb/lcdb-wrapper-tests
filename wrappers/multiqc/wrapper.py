__author__ = "Ryan Dale"
__copyright__ = "Copyright 2016, Ryan Dale"
__email__ = "dalerr@niddk.nih.gov"
__license__ = "MIT"

import os
from snakemake.shell import shell
outdir = os.path.dirname(snakemake.output[0])
if not outdir:
    outdir = '.'

extra = snakemake.params.get('extra', "")
log = snakemake.log_fmt_shell()
shell(
    'multiqc '
    '--quiet '
    '--outdir {outdir} '
    '--force '
    '--filename {snakemake.output} '
    '{extra} '
    '{snakemake.params.analysis_directory} '
    '{log}'
)
