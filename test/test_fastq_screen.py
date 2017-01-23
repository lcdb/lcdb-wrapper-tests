import os
import zipfile
from utils import run, dpath, rm, symlink_in_tempdir

def test_fastqc(sample1_se_fq, tmpdir):
    config = {'database': {
                'ecoli': {
                    'bowtie':
                    },
              },

              'aligner_paths': {
              }
             }
    snakefile = '''
    rule fastq_screen:
        input:
            fastq='sample1_R1.fastq.gz'
        output:
            png='sample1_R1_screen.png',
            txt='sample1_R1_screen.txt'
        params:
            fastq_screen_config={_config},
            subset=100000,
            aligner='bowtie2'
        wrapper: "file://wrapper"'''.format(_config=_config)

    input_data_func=symlink_in_tempdir(
        {
            sample1_se_fq: 'sample1_R1.fastq.gz'
        }
    )

    def check():
        with open('sample1_R1_screen.txt') as fh:
            print(fh.read())

    run(dpath('../wrappers/fastqc'), snakefile, check, input_data_func, tmpdir)
