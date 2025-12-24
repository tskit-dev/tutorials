# tutorials

A set of tutorials for the tskit ecosystem build using
[Jupyter Book](https://jupyterbook.org/), served up at
[https://tskit.dev/tutorials/](https://tskit.dev/tutorials/).

Merges to this repo will trigger a rebuild of the
[tskit.dev web site](https://tskit.dev/) via an
[action](https://github.com/tskit-dev/tskit-site/actions) on the
[tskit-site repository](https://github.com/tskit-dev/tskit-site/):
look there for any deployment issues.

**Under construction**

These are quick notes for developers while the real developers page is
under construction.

# Requirements

Install the Python requirements from `requirements.txt`:
```
$ python -m pip install -r requirements.txt
```

You will also need a working R installation with `reticulate` and `irkernel` installed.
This command should do the trick:
```
$ R -e 'install.packages(c("reticulate", "IRkernel")); IRkernel::installspec()'
```

# Building tutorials

- To add a new tutorial, create a Markdown file and add its name to ``_toc.yml``.
- If you are basing the tutorial on an existing notebook, use
  [jupytext](https://github.com/mwouts/jupytext) to convert the notebook into
  the right format.
- To build locally, run ``make``. The output tells you where to find the
  built `HTML` file(s).
- Pages rendered at https://tskit.dev/tutorials.
- Pages might take a while to be updated after a new tutorial is merged.

If you have an idea for a tutorial, please open an issue to discuss.
