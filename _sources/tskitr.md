---
jupytext:
  formats: md:myst,ipynb
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: R
  language: R
  name: ir
---

(sec_tskit_r)=

# Tskit and R

To interface with `tskit` in R, we can use the [reticulate](https://rstudio.github.io/reticulate/) R package, which lets you call Python functions within an R session. In this short tutorial, we'll go through a couple of examples to show you how to get started. If you haven't done so already, you'll need to install `reticulate` in your R session via `install.packages("reticulate")`. 

We'll begin by simulating a small tree sequence using `msprime`.

```{code-cell}
msprime <- reticulate::import("msprime")

ts <- msprime$sim_ancestry(10, sequence_length = 100, random_seed=42)
ts
```

`reticulate` allows us to access a Python object's attributes or call its methods via the `$` operator. For example, we can access (and assign to a variable) the number of samples in the tree sequence:

```{code-cell}
n <- ts$num_samples
n
```

We can also use `tskit`'s powerful [Statistics](https://tskit.dev/tskit/docs/stable/stats.html) framework to efficiently compute many different summary statistics from a tree sequence. To illustrate this, we'll first add some mutations to our tree sequence with `msprime`'s `sim_mutations` function, and then compute the genetic diversity for each of the tree sequence's sample nodes:

```{code-cell}
ts_mut = msprime$sim_mutations(ts, rate = 1e-2, random_seed = 42)
ts_mut$num_mutations
ts_mut$diversity()
```

As a final example, we can also use the tree sequence `genotype_matrix()` method to return the genotypes of the the tree sequence as a matrix object in R.

```{code-cell}
G = ts_mut$genotype_matrix()
G
```

We can then use R functions directly on the genotype matrix:

```{code-cell}
allele_frequency = rowMeans(G)
allele_frequency
```

It's as simple as that! Be sure to check out the [reticulate](https://rstudio.github.io/reticulate/) documentation, in particular on [Calling Python from R](https://rstudio.github.io/reticulate/articles/calling_python.html), which includes important information on how R data types are converted to their equivalent Python types. 