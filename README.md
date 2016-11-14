# Testing Snakemake wrappers

This is a proof-of-concept for figuring out a good way of testing [Snakemake
wrappers](https://bitbucket.org/snakemake/snakemake-wrappers).

With the [`bioconda`](https://bioconda.github.io/) channel set up, install
the dependencies:

```
conda create -n testenv --file requirements.txt
```

And run tests:

```
source activate testenv
pytest test/tests.py -v
```

## Writing a test

Write new tests in `tests/test_<wrappername>.py`.

A test consists of a Snakefile (typically written as a triple-quoted string
inside each test function), an input data function (typically created via
`utils.symlink_in_tempdir`), and a function that performs arbitrary tests on
the output.

Here's an annotated example.


```python

# (contents of `test/test_fastq_count.py`)

# The input argument may be a test fixture in `conftest.py` that is implicitly
# imported by pytest.
#
# Here we use the `sample1_se_fq` fixture, which will give us the name of the
# tempfile to which a FASTQ file (sample 1's single-end FASTQ from the test
# data repo) has been saved.

# Test utilities we'll be using.
from utils import run, dpath, symlink_in_tempdir

# `run` does most of the work. It creates a tempdir, copys over input data,
# Snakefile, and wrapper, runs the Snakefile, and runs a user-provided test
# function against the output.
#
# `dpath` figures out the path the wrapper
#
# `symlink_in_tempdir` is a decorator function that makes it easy to use pytest
# fixtures that download example data and put it in a place expected by the test.


def test_wrapper(sample1_se_fq):

    # Let's assume we're testing a wrapper that counts the number of fastq reads.

    # Here's the content of our test snakefile with a single rule. The
    # "{wrapper}" placeholder will be filled in automatically with the correct
    # path to the wrapper relative to the tempdir that will be created during the
    # test so we don't have to worry about that here.

    snakefile = '''
                rule test:
                    input: 'sample.fastq.gz'
                    output: 'sample.count'
                    wrapper: {wrapper}
                '''


    # Next we define a function that will be called right before running the
    # Snakefile. Here we use symlink_in_tempdir. We give it a dictionary mapping
    # the temp filename provided by the fixture to its intended location as an
    # input file for the Snakefile.

    input_data_func = symlink_in_tempdir(
        {
            sample1_se_fq: 'sample.fastq.gz'
        }
    )


    # Now we define a function that will test the output. It will be called
    # with the current working directory set to the temporary testing directory,
    # so all paths used in the function are relative to the temp dir.

    def check():
        assert int(open('sample.count')) == 2481



    # Last, we call the `run` function with what we've just defined. Note that
    # paths to wrappers are provided relative to the file from which `dpath` is
    # called.

    run(
        dpath('../wrappers/fastq_count'),
        snakefile=snakefile,
        check=check,
        input_data_func=input_data_func
    )
```

This test, along with any others, can be run using:

```bash
pytest test -v
```
