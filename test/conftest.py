import os
import pytest
import tempfile
import shutil
from snakemake.shell import shell
from snakemake.utils import makedirs
from lcdblib.snakemake import aligners

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
def hisat2_indexes(tmpdir_factory):
    d = str(tmpdir_factory.mktemp('hisat2_indexes'))
    fns = []
    for fn in aligners.hisat2_index_from_prefix('seq/2L'):
        fns.append(_download_file(fn, d))
    return fns


@pytest.fixture(scope='session')
def annotation(tmpdir_factory):
    fn = 'annotation/dm6.gtf'
    d = str(tmpdir_factory.mktemp('annotation'))
    return _download_file(fn, d)
