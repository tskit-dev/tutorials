---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{currentmodule} tskit
```

```{code-cell} ipython3
:tags: [remove-cell]
import msprime
import numpy as np
import tskit

def create_notebook_data():
    pass

# create_notebook_data()  # uncomment to recreate the tree seqs used in this notebook
```

(sec_parallelization)=

# _Parallelization_
% remove underscores in title when tutorial is complete or near-complete

When performing large calculations it's often useful to split the
work over multiple processes or threads. The ``tskit`` API can
be used without issues across multiple processes, and the Python
{mod}`multiprocessing` module often provides a very effective way to
work with many replicate simulations in parallel.

When we wish to work with a single very large dataset, however, threads can
offer better resource usage because of the shared memory space. The Python
{mod}`threading` library gives a very simple interface to lightweight CPU
threads and allows us to perform several CPU intensive tasks in parallel. The
``tskit`` API is designed to allow multiple threads to work in parallel when
CPU intensive tasks are being undertaken.

:::{note}
In the CPython implementation the 
[Global Interpreter Lock](https://wiki.python.org/moin/GlobalInterpreterLock) ensures that
only one thread executes Python bytecode at one time. This means that
Python code does not parallelise well across threads, but avoids a large
number of nasty pitfalls associated with multiple threads updating
data structures in parallel. Native C extensions like ``numpy`` and ``tskit``
release the GIL while expensive tasks are being performed, therefore
allowing these calculations to proceed in parallel.
:::


:::{todo}
This tutorial previously used code with an old interface, and hence has been removed.
We must recreate an example of parallel processing, giving examples of both
threads and processes (but see
[this stackoverflow post](https://stackoverflow.com/questions/47313732/jupyter-notebook-never-finishes-processing-using-multiprocessing-python-3)
for why it may be difficult to get {mod}`multiprocessing` working in this notebook).
A reasonable example might be to calculate many pairwise statistics between sample sets
in parallel.

We should also show how, for large tree sequences that it is better to pass the filenames
to each subprocess, and load the tree sequence, rather than transferring the entire
tree sequence (via pickle) to the subprocesses.
:::

