import nbformat
import nbconvert
from nbconvert.preprocessors import ExecutePreprocessor
import sys

nbfile = sys.argv[1]
nboutfile = sys.argv[2]
with open(nbfile) as f:
    nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=600, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': '.'}})
    with open(nboutfile, 'wt') as f:
        nbformat.write(nb, f)
