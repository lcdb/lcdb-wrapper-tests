import os
import pytest
import tempfile
import shutil
from snakemake.shell import shell
from snakemake.utils import makedirs
from lcdblib.snakemake import aligners
from utils import run, dpath, symlink_in_tempdir

# test data url
URL = 'https://github.com/lcdb/lcdb-test-data/blob/master/data/{}?raw=true'

def _download_file(fn, d):
    """
    Intended to be called from a pytest.fixture function below.

    `fn` is a path to a file that is used to fill in `URL`. `d` is a tempdir
    likely created by the calling function to which the file will be
    downloaded.

    The path to the downloaded file is returned.
    """
    url = URL.format(fn)
    dest = os.path.join(d, fn)
    makedirs(os.path.dirname(dest))
    basename = os.path.basename(fn)
    shell('wget -q -O- {url} > {dest}')
    return dest


@pytest.fixture(scope='session')
def sample1_se_fq(tmpdir_factory):
    d = str(tmpdir_factory.mktemp('sample1_se_fq'))
    fn = 'samples/sample1/sample1_R1.fastq.gz'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_pe_fq(tmpdir_factory):
    pair = []
    d = str(tmpdir_factory.mktemp('sample1_pe_fq'))
    for fn in [
        'samples/sample1/sample1_R1.fastq.gz',
        'samples/sample1/sample1_R2.fastq.gz'
    ]:
        pair.append(_download_file(fn, d))
    return pair


@pytest.fixture(scope='session')
def dm6_fa(tmpdir_factory):
    fn = 'seq/2L.fa'
    d = str(tmpdir_factory.mktemp('dm6_fa'))
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_se_bam(tmpdir_factory):
    d = str(tmpdir_factory.mktemp('sample1_se_bam'))
    fn = 'samples/sample1/sample1.single.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_pe_bam(tmpdir_factory):
    d = str(tmpdir_factory.mktemp('sample1_pe_bam'))
    fn = 'samples/sample1/sample1.paired.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def hisat2_indexes(dm6_fa, tmpdir_factory):
    d = str(tmpdir_factory.mktemp('hisat2_indexes'))
    snakefile = '''
    rule hisat2:
        input: fasta='2L.fa'
        output: index=['2L.1.ht2', '2L.2.ht2']
        wrapper: 'file://wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            dm6_fa: '2L.fa'
        }
    )
    def check():
        pass
    run(dpath('../wrappers/hisat2/build'), snakefile, check, input_data_func, d)
    return aligners.hisat2_index_from_prefix(os.path.join(d, '2L'))


@pytest.fixture(scope='session')
def bowtie2_indexes(dm6_fa, tmpdir_factory):
    d = str(tmpdir_factory.mktemp('bowtie2_indexes'))
    snakefile = '''
    rule bowtie2:
        input: fasta='2L.fa'
        output: index=['2L.1.bt2', '2L.2.bt2']
        wrapper: 'file://wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            dm6_fa: '2L.fa'
        }
    )
    def check():
        pass
    run(dpath('../wrappers/bowtie2/build'), snakefile, check, input_data_func, d)
    return aligners.bowtie2_index_from_prefix(os.path.join(d, '2L'))

@pytest.fixture(scope='session')
def annotation(tmpdir_factory):
    fn = 'annotation/dm6.gtf'
    d = str(tmpdir_factory.mktemp('annotation'))
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def annotation_refflat(tmpdir_factory):
    fn = 'annotation/dm6.refflat'
    d = str(tmpdir_factory.mktemp('annotation_refflat'))
    return _download_file(fn, d)

@pytest.fixture(scope='session')
def fastqc(sample1_se_fq, tmpdir_factory):
    snakefile = '''
    rule fastqc:
        input:
            fastq='sample1_R1.fastq.gz'
        output:
            html='sample1_R1_fastqc.html',
            zip='sample1_R1_fastqc.zip'
        wrapper: "file://wrapper"'''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_fq: 'sample1_R1.fastq.gz'
        }
    )
    tmpdir = str(tmpdir_factory.mktemp('fastqc_fixture'))
    run(dpath('../wrappers/fastqc'), snakefile, None, input_data_func, tmpdir)
    return os.path.join(tmpdir, 'sample1_R1_fastqc.zip')
