import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir

def test_fastqc(sample1_se_fq, tmpdir):
    snakefile = '''
    rule fastqc:
        input:
            fastq='sample1_R1.fastq.gz'
        output:
            html='results/sample1_R1.html',
            zip='sample1_R1.zip'
        wrapper: "file://wrapper"'''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        assert '<html>' in open('results/sample1_R1.html').readline()

    run(dpath('../wrappers/fastqc'), snakefile, check, input_data_func, tmpdir)
