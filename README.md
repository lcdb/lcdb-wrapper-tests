# Testing Snakemake wrappers

This is a proof-of-concept for figuring out a good way of testing [Snakemake
wrappers](https://bitbucket.org/snakemake/snakemake-wrappers).

With the [`bioconda`](https://bioconda.github.io/) channel set up, install
the dependencies:

```
conda create -n testenv snakemake pytest
```

And run tests:

```
source activate testenv
pytest test/tests.py -v
```

## Writing a test

In order to run a test:

- The wrapper directory must contain a `test` directory.
- Inside that must be one or more arbitrarily-named test-case directories.
- Each test-case directory must have at least a Snakefile and an
expected-results subdirectory. If the test uses any input files, they must be
in the test-case directory.

For example, here is a wrapper that has one test-case called "test1":

```
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
```
For each test-case, a tempdir will be created and the wrapper contents will
be copied there (environment.yml, README.md, and wrapper.py). The test data
will be copied to that same directory. That is, for test-case "test1"
above, the test will be run as `snakemake --workdir TMPDIR` and the
directory structure looks like this:
```
	TMPDIR/
	├── environment.yml
	├── README.md
	├── Snakefile
	├── testfile
	└── wrapper.py
```
Upon running the test, in this example we expect the `output` file to be
created, so a successful test will have the following directory structure:
```

	TMPDIR/
	├── environment.yml
	├── output
	├── README.md
	├── Snakefile
	├── stats.txt
	├── testfile
	└── wrapper.py
```

This means that in the test Snakefile, the wrapper can refer to the current
working directory as "file://". A simple test Snakefile might look like
this:

```
	rule test:
		input: 'testfile'
		output: 'output'
		wrapper: 'file://'
```
