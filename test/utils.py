"""
Stripped-down version of Snakemake's test framework.
"""

import sys
import os
import subprocess as sp
import tempfile
import hashlib
import urllib
import shutil
import shlex

import pytest
from snakemake import snakemake


SCRIPTPATH = shutil.which('snakemake')


def dpath(path):
    "path relative to this file"
    return os.path.realpath(os.path.join(os.path.dirname(__file__), path))


def md5sum(filename):
    data = open(filename, 'rb').read()
    return hashlib.md5(data).hexdigest()


def run(path, test_case, check_md5=True, input_data_func=None, output_data_func=None,
        **params):
    """
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

    """
    # store any tempdirs here for later deletion
    to_clean_up = []

    test_path = os.path.join(path, 'test', test_case)

    # Call the output_data_func if provided, otherwise we should find the
    # expected-results dir
    if output_data_func is not None:
        results_dir = tempdir.mkdtemp(prefix='.expected-results', dir=os.path.abspath('.'))
        to_clean_up.append(results_dir)
        output_data_func(results_dir)
    else:
        results_dir = os.path.join(test_path, 'expected-results')
    assert (
        os.path.exists(results_dir) and os.path.isdir(results_dir)
    ), '{} does not exist'.format(results_dir)

    # Should always find a Snakefile
    snakefile = os.path.join(test_path, 'Snakefile')
    assert os.path.exists(snakefile), 'expected {}'.format(snakefile)

    tmpdir = tempfile.mkdtemp(prefix='.test', dir=os.path.abspath('.'))
    try:
        # Copy over the test files
        cmds = (
            'find {} -maxdepth 1 -type f -print0 | xargs -0 cp -t {}'
            .format(shlex.quote(test_path), shlex.quote(tmpdir))
        )
        sp.call(cmds, shell=True)

        # Copy over the wrapper files
        cmds = (
            'find {} -maxdepth 1 -type f -print0 | xargs -0 cp -t {}'
            .format(shlex.quote(path), shlex.quote(tmpdir))
        )
        sp.call(cmds, shell=True)

        # Create the input data if needed
        if input_data_func is not None:
            input_data_func(tmpdir)

        success = snakemake(snakefile, workdir=tmpdir, stats='stats.txt',
                            snakemakepath=SCRIPTPATH, config={}, **params)
        assert success, 'expected successful execution'

        # Check that for everything in the results_dir we got a matching file
        # created in the tmpdir, comparing md5sums if specified.
        for resultfile in os.listdir(results_dir):

            # skip directories or known-to-be-uninteresting files
            if (
                resultfile == '.gitignore' or
                not os.path.isfile(
                    os.path.join(results_dir, resultfile)
                )
            ):
                continue

            targetfile = os.path.join(tmpdir, resultfile)
            expectedfile = os.path.join(results_dir, resultfile)
            assert os.path.exists(
                targetfile), 'expected file "{}" not created'.format(
                    resultfile)
            if check_md5:
                assert md5sum(targetfile) == md5sum(
                    expectedfile), 'md5sum mismatch for "{}"'.format(
                        resultfile)
    finally:
        for t in to_clean_up:
            shutil.rmtree(t)
