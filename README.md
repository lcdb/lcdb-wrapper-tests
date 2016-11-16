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

## How the fixtures work

`py.test` *fixtures* are functions. They are used for things like setup and
teardown that needs to be done once, or for preparing data just once that can
be used across many tests. In this repo, the fixtures are primarily used for
downloading data from the testing data repo to a test-specific temporary
directory.

Fixtures are stored in a `conftest.py` file (see `tests/conftest.py`). Any
discovered pytest tests in any module (functions starting with `test_`) can use
any of these fixtures as an argument, and will get the return value of that
fixture at runtime.

pytest also has some built-in fixtures. Here we use the `tmpdir` fixture and
the `tmpdir_factory` fixture. The nice thing about these is that pytest
manages them on disk and only keeps the last handful of tempdirs. This makes it
straightforward to find the test data on disk and poke around when
troubleshooting.

### Specific example from the test suite

Let's take a concrete example. In the `tests/conftest.py` file, there's
a function called `sample1_se_fq`. Its job is to download a single FASTQ file
from the test data repo. It is marked with a `pytest.fixture` decorator, and is
set to run only once per session. Fixtures can use other fixtures, and
`sample1_se_fq()` uses the built-in `tmpdir_factory` fixture which is used for
once-per-session fixtures like we have here.

Once per test session (that is, a single invocation of the `pytest` command),
the `sample1_se_fq()` fixture will download the single-end FASTQ file to
`/tmp/pytest-of-$USER/pytest-$N/sample1_se_fq0/` (filling in the current user
and the current test number N). The fixture stores the path it downloaded and
can provide it to other fixtures or tests.

Over in `tests/test_cutadapt.py`, we have tests for the cutadapt wrapper. For
example the `test_cutadapt_simple` test has two fixtures: `sample1_se_fq` and
`tmpdir`. So it will get the path to that downloaded FASTQ as its first
argument, and will get the built-in pytest fixture as its second argument. You
can read the code for details but basically it copies the Snakefile and
configured input data (via the `symlink_in_tempdir` function) to the directory
`/tmp/pytest-of-$USER/pytest-$N/test_cutadapt_simple0/`.

**Importantly, this means that if a test fails and you need to troubleshoot,
you can go to `/tmp/pytest-of-$USER/pytest-$N/$TEST_NAME` directory and you'll
see the Snakefile, input data, and any created output files.**
