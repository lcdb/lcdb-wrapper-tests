import os
import pytest
import gzip
import tempfile
import shutil
from utils import run, dpath, rm, symlink_in_tempdir

def test_cutadapt_simple(sample1_se_fq):
    snakefile = '''
                rule cutadapt:
                    input:
                        fastq='sample1_R1.fastq.gz'
                    output:
                        fastq='sample1_R1.trim.fastq.gz'
                    params: extra='-a AAA'
                    wrapper: {wrapper}
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        a = sum(1 for _ in gzip.open('sample1_R1.fastq.gz'))
        b = sum(1 for _ in gzip.open('sample1_R1.trim.fastq.gz'))
        assert a == b == 9924

        assert os.path.getsize('sample1_R1.fastq.gz') != os.path.getsize('sample1_R1.trim.fastq.gz')

    run(dpath('../wrappers/cutadapt'), snakefile, check, input_data_func)


def test_cutadapt_simple_with_log(sample1_se_fq):
    snakefile = '''
                rule cutadapt:
                    input:
                        fastq='sample1_R1.fastq.gz'
                    output:
                        fastq='sample1_R1.trim.fastq.gz'
                    params: extra='-a AAA'
                    log: 'sample1.cutadapt.log'
                    wrapper: {wrapper}
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        a = sum(1 for _ in gzip.open('sample1_R1.fastq.gz'))
        b = sum(1 for _ in gzip.open('sample1_R1.trim.fastq.gz'))
        assert a == b == 9924
        assert 'This is cutadapt' in open('sample1.cutadapt.log').readline()

        assert os.path.getsize('sample1_R1.fastq.gz') != os.path.getsize('sample1_R1.trim.fastq.gz')

    run(dpath('../wrappers/cutadapt'), snakefile, check, input_data_func)

