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

(sec_incremental)=
# _Incremental algorithms_

Much of the [efficiency](sec_what_is_analysis)
of the tskit approach comes from the use of incremental algorithms.
By considering only the difference between adjacent trees,
incremental algorithms avoid having to perform the same
calculation multiple times on different trees.

This tutorial will explain the philosophy behind incremental algorithms,
and provide examples of how to create your own (e.g. using the
{meth}`TreeSequence.edge_diffs` method).

:::{todo}
Create content. See [issue 233](https://github.com/tskit-dev/tutorials/issues/233)
:::