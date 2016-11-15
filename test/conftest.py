import os
import pytest
import tempfile
import shutil
from snakemake.shell import shell

@pytest.fixture(scope='session')
def sample1_se_fq(tmpdir_factory):
    fn = 'samples/sample1/sample1_R1.fastq.gz'
    d = str(tmpdir_factory.mktemp('sample1_se_fq'))
    url = 'https://github.com/lcdb/lcdb-test-data/blob/add-data/data/{}?raw=true'.format(fn)
    basename = os.path.basename(fn)
    shell('wget -O- {url} > {d}/{basename}')
    return os.path.join(d, basename)


@pytest.fixture(scope='session')
def sample1_pe_fq(tmpdir_factory):
    pair = []
    d = str(tmpdir_factory.mktemp('sample1_pe_fq'))
    for fn in [
        'samples/sample1/sample1_R1.fastq.gz',
        'samples/sample1/sample1_R2.fastq.gz'
    ]:
        url = 'https://github.com/lcdb/lcdb-test-data/blob/add-data/data/{}?raw=true'.format(fn)
        basename = os.path.basename(fn)
        shell('wget -O- {url} > {d}/{basename}')
        pair.append(os.path.join(d, basename))
    return pair
