import os
import pytest
import tempfile
import shutil
from snakemake.shell import shell
from lcdblib.snakemake import aligners

@pytest.fixture(scope='session')
def sample1_se_fq(tmpdir_factory):
    fn = 'samples/sample1/sample1_R1.fastq.gz'
    d = str(tmpdir_factory.mktemp('sample1_se_fq'))
    url = 'https://github.com/lcdb/lcdb-test-data/blob/master/data/{}?raw=true'.format(fn)
    basename = os.path.basename(fn)
    shell('wget -q -O- {url} > {d}/{basename}')
    return os.path.join(d, basename)


@pytest.fixture(scope='session')
def sample1_pe_fq(tmpdir_factory):
    pair = []
    d = str(tmpdir_factory.mktemp('sample1_pe_fq'))
    for fn in [
        'samples/sample1/sample1_R1.fastq.gz',
        'samples/sample1/sample1_R2.fastq.gz'
    ]:
        url = 'https://github.com/lcdb/lcdb-test-data/blob/master/data/{}?raw=true'.format(fn)
        basename = os.path.basename(fn)
        shell('wget -q -O- {url} > {d}/{basename}')
        pair.append(os.path.join(d, basename))
    return pair


@pytest.fixture(scope='session')
def dm6_fa(tmpdir_factory):
    fn = 'seq/2L.fa'
    d = str(tmpdir_factory.mktemp('dm6_fa'))
    url = 'https://github.com/lcdb/lcdb-test-data/blob/master/data/{}?raw=true'.format(fn)
    basename = os.path.basename(fn)
    shell('wget -q -O- {url} > {d}/{basename}')
    return os.path.join(d, basename)


@pytest.fixture(scope='session')
def hisat2_indexes(tmpdir_factory):
    d = str(tmpdir_factory.mktemp('hisat2_indexes'))
    fns = []
    for fn in aligners.hisat2_index_from_prefix('seq/2L'):
        url = 'https://github.com/lcdb/lcdb-test-data/blob/master/data/{}?raw=true'.format(fn)
        basename = os.path.basename(fn)
        shell('wget -q -O- {url} > {d}/{basename}')
        fns.append(os.path.join(d, basename))
    return fns
