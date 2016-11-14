import os
import pytest
import tempfile
import shutil
from snakemake.shell import shell

@pytest.fixture(scope='module')
def sample1_se_fq():
    fn = 'samples/sample1/sample1_R1.fastq.gz'
    d = tempfile.mkdtemp()
    url = 'https://github.com/lcdb/lcdb-test-data/blob/add-data/data/{}?raw=true'.format(fn)
    basename = os.path.basename(fn)
    shell('wget -O- {url} > {d}/{basename}')
    yield os.path.join(d, basename)
    shutil.rmtree(d)
