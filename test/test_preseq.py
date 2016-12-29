import os
import gzip
from utils import run, dpath, symlink_in_tempdir
from textwrap import dedent

def test_infer_experiment(sample1_se_sort_bam, tmpdir):
    snakefile = '''
                rule observed_complexity:
                    input:
                        bam='sample1_R1.sort.bam'
                    output:
                        table='sample1_R1.sort.observed_complexity.txt'
                    wrapper: "file://wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_sort_bam: 'sample1_R1.sort.bam'
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        expected = ('total_reads\tdistinct_reads\n'
                    '0\t0')

        with open('sample1_R1.sort.observed_complexity.txt', 'r') as handle:
            results = handle.read().strip()

        assert  results == expected

    run(dpath('../wrappers/preseq/observed_complexity'), snakefile, check,
              input_data_func, tmpdir, use_conda=True)

