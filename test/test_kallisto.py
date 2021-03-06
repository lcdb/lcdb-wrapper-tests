import pytest
import json
import pysam
from snakemake.shell import shell
from lcdblib.snakemake import aligners
from utils import run, dpath, rm, symlink_in_tempdir


def test_kallisto_index(transcriptome, tmpdir):
    snakefile = '''
    rule kallisto:
        input: fasta='transcriptome.fa'
        output: index='out/transcriptome.idx'
        log: 'log'
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            transcriptome: 'transcriptome.fa',
        }
    )

    def check():
        log = open('log').read()
        assert '[build] target deBruijn graph'

    run(
        dpath('../wrappers/kallisto/index'),
        snakefile, check, input_data_func, tmpdir)


def test_kallisto_quant(tmpdir, sample1_se_fq, kallisto_index):
    snakefile = '''
    rule kallisto_quant:
        input:
             fastq='sample1.fq.gz',
             index='out/transcriptome.idx'

        params: extra='--single --fragment-length=200 --sd=20'
        output:
            h5='quant/abundance.h5',
            tsv='quant/abundance.tsv',
            json='quant/run_info.json',
        wrapper: 'file:wrapper'
    '''
    input_data_func = symlink_in_tempdir(
        {
            sample1_se_fq: 'sample1.fq.gz',
            kallisto_index: 'out/transcriptome.idx',
        }
    )

    def check():
        assert sum(1 for _ in open('quant/abundance.tsv')) == 35
        assert open('quant/abundance.tsv').readline() == (
                'target_id\tlength\teff_length\test_counts\ttpm\n')
        keys = ['call', 'index_version', 'n_bootstraps', 'n_processed', 'n_targets', 'start_time']
        d = json.load(open('quant/run_info.json'))
        for k in keys:
            assert k in d


    run(
        dpath('../wrappers/kallisto/quant'),
        snakefile, check, input_data_func, tmpdir)
