from snakemake.shell import shell
from utils import run, dpath, rm, symlink_in_tempdir


def test_hisat2_build(dm6_fa, tmpdir):
    snakefile = '''
                rule hisat2_build:
                    input:
                        fasta='2L.fa'
                    output:
                        index=expand('data/assembly/assembly.{n}.ht2', n=range(1,9))
                    log: 'hisat.log'
                    wrapper: "file://wrapper"

                '''
    input_data_func=symlink_in_tempdir(
        {
            dm6_fa: '2L.fa'
        }
    )

    def check():
        assert 'Total time for call to driver' in open('hisat.log').readlines()[-1]
        assert list(shell('hisat2-inspect data/assembly/assembly -n', iterable=True)) == ['2L']

    run(dpath('../wrappers/hisat2/build'), snakefile, check, input_data_func, tmpdir)
