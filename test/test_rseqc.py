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


def test_geneBody_cov(sample1_se_sort_bam, sample1_se_sort_bam_bai, annotation_bed12, tmpdir):
    snakefile = '''
                rule geneBody_coverage:
                    input:
                        bam='sample1_R1.sort.bam',
                        bai='sample1_R1.sort.bam.bai',
                        bed='dm6.bed12'
                    output: txt='sample1_R1.geneBodyCoverage.txt',
                            r='sample1_R1.geneBodyCoverage.r',
                            pdf='sample1_R1.geneBodyCoverage.pdf',
                    wrapper: "file://wrapper"
                '''
    input_data_func=symlink_in_tempdir(
        {
            sample1_se_sort_bam: 'sample1_R1.sort.bam',
            sample1_se_sort_bam_bai: 'sample1_R1.sort.bam.bai',
            annotation_bed12: 'dm6.bed12'
        }
    )

    def check():
        """
        check for line lengths and that they are at least different sized
        """

        # R code
        with open('sample1_R1.geneBodyCoverage.r', 'r') as handle:
            result = handle.readline().split(' ')[0]

        assert  result == 'sample1_R1.sort'

        # text
        with open('sample1_R1.geneBodyCoverage.txt', 'r') as handle:
            result = handle.readlines()[1].split('\t')[0]

        assert  result == 'sample1_R1.sort'

        # PDF
        assert os.path.exists('sample1_R1.geneBodyCoverage.pdf')

    run(dpath('../wrappers/rseqc/geneBody_covearge'), snakefile, check, input_data_func, tmpdir, use_conda=True)

