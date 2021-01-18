import nbformat
import nbconvert
from nbconvert.preprocessors import ExecutePreprocessor
import argparse
import sys


def make_parser():
    parser = argparse.ArgumentParser(
        description="Options for converting notebooks")

    parser.add_argument('--timeout', type=int, default=600,
                        help='Execution timeout (seconds)')
    return parser


nbfile = sys.argv[1]
nboutfile = sys.argv[2]
parser=make_parser()
args=parser.parse_args(sys.argv[3:])
with open(nbfile) as f:
    nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=args.timeout, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': '.'}})
    with open(nboutfile, 'wt') as f:
        nbformat.write(nb, f)
