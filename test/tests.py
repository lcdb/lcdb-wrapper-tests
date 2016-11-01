import os
import pytest
from utils import run, dpath

def test_ls():
    run(dpath('../wrappers/ls'), 'test1')


def ls_input_data(path):
    with open(os.path.join(path, 'test2_infile'), 'w') as fout:
        fout.write('contents')

@pytest.mark.xfail
def test_ls2_fails():
    """
    This test doesn't provide its own input data and relies on ls_input_data,
    so not providing that function should fail
    """
    run(dpath('../wrappers/ls'), 'test2')

def test_ls2():
    run(dpath('../wrappers/ls'), 'test2', input_data_func=ls_input_data)

