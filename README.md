# Testing Snakemake wrappers

This is a proof-of-concept for figuring out a good way of testing [Snakemake
wrappers](https://bitbucket.org/snakemake/snakemake-wrappers).

With the [`bioconda`](https://bioconda.github.io/) channel set up, install
dependencies:

```
conda create -n testenv snakemake pytest
```

And run tests:

```
source activate testenv
pytest test/tests.py -v
```

Most of the work happens in `test.utils.run`, so the docstring of that function is copied here:

```
    Parameters
    ----------

    path : str
        Path to a wrapper directory.

    test_case : str
        Name of the test-case, see below.

    check_md5 : bool
        If True, then ensure the md5sums of generated output files and expected
        output files match.

    input_data_func : None | callable
        If not None, then this callable object will be called with
        a single argument corresponding to the temp directory. It will be
        called after the wrapper and test-case contents have been copied to the
        temp dir, but before the test is run. It is expected to create any data
        required in whatever directory structure is required.

    output_data_func : None | callable
        If not None, then this callable object will be called with a single
        argument corresponding to temporary directory, which is separate from
        the temp dir created for running the test. This new temp dir works as
        a replacement for the expected-results directory.

    params :
        Extra kwargs passed to the snakemake function.

    In order to run a test, the wrapper directory must contain a "test"
    directory. Inside that must be one or more arbitrarily-named test-case
    directories. Each test-case directory must have at least a Snakefile and an
    expected-results subdirectory. If the test uses any input files, they must
    be in the test-case directory.

    For example, here is a wrapper that has one test-case called "test1"::

        WRAPPER_DIRECTORY
        ├── environment.yml
        ├── README.md
        ├── test
        │   └── test1
        │       ├── expected-results
        │       │   └── output
        │       ├── Snakefile
        │       └── testfile
        └── wrapper.py

    For each test-case, a tempdir will be created and the wrapper contents will
    be copied there (environment.yml, README.md, and wrapper.py). The test data
    will be copied to that same directory. That is, for test-case "test1"
    above, the test will be run as `snakemake --workdir TMPDIR` and the
    directory structure looks like this::

        TMPDIR/
        ├── environment.yml
        ├── README.md
        ├── Snakefile
        ├── testfile
        └── wrapper.py

    Upon running the test, in this example we expect the `output` file to be
    created, so a successful test will have the following directory structure::

        TMPDIR/
        ├── environment.yml
        ├── output
        ├── README.md
        ├── Snakefile
        ├── stats.txt
        ├── testfile
        └── wrapper.py

    This means that in the test Snakefile, the wrapper can refer to the current
    working directory as "file://". A simple test Snakefile might look like
    this::

        rule test:
            input: 'testfile'
            output: 'output'
            wrapper: 'file://'
```
