import os
import gzip
from utils import run, dpath, rm, symlink_in_tempdir

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
        with open('sample1_R1.infer_experiment.txt', 'r') as handle:
            line = handle.readline()
            assert line == 'This is SingleEnd Data'

    run(dpath('../wrappers/rseqc/infer_experiment'), snakefile, check, input_data_func, tmpdir, use_conda=True)
