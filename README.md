[![Build Status](https://travis-ci.org/lcdb/lcdb-wrapper-tests.svg?branch=master)](https://travis-ci.org/lcdb/lcdb-wrapper-tests)

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

See `test/test_demo.py` for an annotated example.
