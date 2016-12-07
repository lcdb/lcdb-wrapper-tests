import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir
from textwrap import dedent

def test_infer_experiment(sample1_se_bam, annotation_bed12, tmpdir):
    snakefile = '''
                rule infer_experiment:
                    input:
                        bam='sample1_R1.bam',
                        bed='dm6.bed12'
                    output: 'sample1_R1.infer_experiment.txt'
                    wrapper: "file://wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_bam: 'sample1_R1.bam',
            annotation_bed12: 'dm6.bed12'
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """
        expected = dedent("""\
                This is SingleEnd Data
                Fraction of reads failed to determine: 0.1166
                Fraction of reads explained by "++,--": 0.8820
                Fraction of reads explained by "+-,-+": 0.0014""")

        with open('sample1_R1.infer_experiment.txt', 'r') as handle:
            results = handle.read().strip()

        assert  results == expected

    run(dpath('../wrappers/rseqc/infer_experiment'), snakefile, check, input_data_func, tmpdir, use_conda=True)
