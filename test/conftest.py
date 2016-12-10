import os
import pytest
import tempfile
import shutil
import inspect
from snakemake.shell import shell
from snakemake.utils import makedirs
from lcdblib.snakemake import aligners
from utils import run, dpath, symlink_in_tempdir

# test data url
URL = 'https://github.com/lcdb/lcdb-test-data/blob/master/data/{}?raw=true'


def tmpdir_for_func(factory):
    caller = inspect.stack()[1][3]
    return str(factory.mktemp(caller))


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
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'samples/sample1/sample1_R1.fastq.gz'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_pe_fq(tmpdir_factory):
    pair = []
    d = tmpdir_for_func(tmpdir_factory)
    for fn in [
        'samples/sample1/sample1_R1.fastq.gz',
        'samples/sample1/sample1_R2.fastq.gz'
    ]:
        pair.append(_download_file(fn, d))
    return pair


@pytest.fixture(scope='session')
def dm6_fa(tmpdir_factory):
    fn = 'seq/2L.fa'
    d = tmpdir_for_func(tmpdir_factory)
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_se_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'samples/sample1/sample1.single.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_se_sort_bam(sample1_se_bam):
    sample1_se_sort_bam = sample1_se_bam.replace('.bam', '.sort.bam')
    shell(
            "samtools sort "
            "-o {sample1_se_sort_bam} "
            "-O BAM "
            "{sample1_se_bam} "
            )
    return sample1_se_sort_bam


@pytest.fixture(scope='session')
def sample1_se_sort_bam_bai(sample1_se_sort_bam):
    shell(
            "samtools index "
            "{sample1_se_sort_bam}"
            )
    return sample1_se_sort_bam + '.bai'


@pytest.fixture(scope='session')
def sample1_pe_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
    fn = 'samples/sample1/sample1.paired.bam'
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def sample1_pe_hisat2_bam(tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)


@pytest.fixture(scope='session')
def hisat2_indexes(dm6_fa, tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
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

    run(
        dpath('../wrappers/hisat2/build'),
        snakefile, check, input_data_func, d)
    return aligners.hisat2_index_from_prefix(os.path.join(d, '2L'))


@pytest.fixture(scope='session')
def bowtie2_indexes(dm6_fa, tmpdir_factory):
    d = tmpdir_for_func(tmpdir_factory)
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

    run(
        dpath('../wrappers/bowtie2/build'),
        snakefile, check, input_data_func, d)
    return aligners.bowtie2_index_from_prefix(os.path.join(d, '2L'))


@pytest.fixture(scope='session')
def annotation(tmpdir_factory):
    fn = 'annotation/dm6.gtf'
    d = tmpdir_for_func(tmpdir_factory)
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def annotation_refflat(tmpdir_factory):
    fn = 'annotation/dm6.refflat'
    d = tmpdir_for_func(tmpdir_factory)
    return _download_file(fn, d)


@pytest.fixture(scope='session')
def annotation_db(annotation):
    import gffutils
    gffutils.create_db(
        data=annotation, dbfn=annotation + '.db',
        merge_strategy='merge',
        id_spec={'transcript': ['transcript_id', 'transcript_symbol'],
                 'gene': ['gene_id', 'gene_symbol']},
        gtf_transcript_key='transcript_id',
        gtf_gene_key='gene_id')

    return annotation + '.db'


@pytest.fixture(scope='session')
def annotation_bed12(annotation_db):
    import gffutils
    db = gffutils.FeatureDB(annotation_db)

    bed12 = '.'.join(annotation_db.strip().split('.')[:-2]) + '.bed12'

    with open(bed12, 'w') as handle:
        for t in db.features_of_type('transcript'):
            handle.write(db.bed12(t, name_field='transcript_id') + '\n')

    return bed12


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
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_fq: 'sample1_R1.fastq.gz'
        }
    )
    tmpdir = str(tmpdir_factory.mktemp('fastqc_fixture'))
    run(dpath('../wrappers/fastqc'), snakefile, None, input_data_func, tmpdir)
    return os.path.join(tmpdir, 'sample1_R1_fastqc.zip')
