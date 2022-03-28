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

(sec_intro)=

# Welcome!

This site contains a number of tutorials to develop your understanding of
[succinct tree sequences](https://tskit.dev/learn.html#what) and software programs,
such as [msprime](https://tskit.dev/msprime), that use them.

If you are new to the world of tree sequences, we suggest you start with the
first tutorial: {ref}`sec_what_is`

:::{note}
Tutorials are under constant development. Those that are still a work in progress and
not yet ready for use are shown in _italics_ in the list of tutorials.

We very much welcome help developing existing tutorials or writing new ones. Please open
or contribute to a [GitHub issue](https://github.com/tskit-dev/tutorials/issues) if you
would like to help out.
:::

## Other sources of help

In addition to these tutorials, our [Learn page](https://tskit.dev/learn.html) lists
selected videos and publications to help you learn about tree sequences. 

We aim to be a friendly, welcoming open source community.
Questions and discussion about using {program}`tskit`, the tree sequence toolkit
should be directed to the
[GitHub discussion forum](https://github.com/tskit-dev/tskit/discussions), and there are
similar forums for other software in the tree sequence [development community](https://github.com/tskit-dev),
such as for [msprime](https://github.com/tskit-dev/msprime/discussions) and
[tsinfer](https://github.com/tskit-dev/tsinfer/discussions).


(sec_intro_running)=

## Running tutorial code

It is possible to run the tutorial code on your own computer, if you wish.
This will allow you to experiment with the examples provided.
The recommended way to do this is from within a
[Jupyter notebook](https://jupyter.org). As well as installing Jupyter, you will also
need to install the various Python libraries, most importantly
``tskit``, ``msprime``, ``numpy``, and ``matplotlib``. These and other packages are
listed  in the [requirements.txt](https://tskit.dev/tutorials/requirements.txt)
file; a shortcut to installing the necessary software is therefore:

```
python3 -m pip install -r https://tskit.dev/tutorials/requirements.txt
```

In addition, to run the {ref}`R tutorial<sec_tskit_r>` you will need to install the R
[reticulate](https://rstudio.github.io/reticulate/) library, and if running it in a Jupyter
notebook, the [IRkernel](https://irkernel.github.io) library. This can be done by
running the following command within R:

```
install.packages(c("reticulate", "IRkernel")); IRkernel::installspec()
```

(sec_intro_downloading_datafiles)=

### Downloading tutorial datafiles

Many of the tutorials use pre-existing tree sequences stored in the
[``data``](https://github.com/tskit-dev/tutorials/tree/main/data) directory.
These can be downloaded individually from that link, or you can
download them all at once by running the script stored in
[https://tskit.dev/tutorials/examples/download.py](https://tskit.dev/tutorials/examples/download.py).
If you are running the code in the tutorials from within a Jupyter notebook
then you can simply load this code into a new cell by using the
[%load cell magic](https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-load).
Just run the following in a Jupyter code cell:

```
%load https://tskit.dev/tutorials/examples/download.py
```

Running the resulting Python code should download the data files, then print out
``finished downloading`` when all files are downloaded. You should then be able
to successfully run code such as the following:

```{code-cell} ipython3
import tskit
ts = tskit.load("data/basics.trees")
print(f"The file 'data/basics.trees' exists, and contains {ts.num_trees} trees")
```
