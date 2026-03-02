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

```{currentmodule} tskit
```

(sec_tskit_r)=

# Tskit and R

To interface with `tskit` in R, we can use the [reticulate](https://rstudio.github.io/reticulate)
package, which lets you call `tskit`'s Python API from an R session. In this
tutorial, we walk through a few simple examples to help you get started.

First, install `reticulate` in R with `install.packages("reticulate")`. Then,
install the required Python packages with
`reticulate::py_require(c("tskit", "msprime"))`.

By default, this uses `reticulate`'s isolated Python virtual environment, but
you can also use another Python environment as described at
<https://rstudio.github.io/reticulate/>.

```
install.packages("reticulate")
reticulate::py_require(c("tskit", "msprime"))
```

::::{margin}
:::{note}
The [slendr](https://www.slendr.net) R package uses `reticulate` extensively to work with tree sequences in R and also provides R wrappers for some `tskit` functions.
:::
::::

::::{margin}
:::{note}
The [RcppTskit](https://cran.r-project.org/package=RcppTskit) R package provides
access to `tskit`'s C API for advanced use cases and also shows how to interface
with `reticulate` Python.
:::
::::

We'll begin by simulating a small tree sequence
by calling `msprime$sim_ancestry()` from R.

```{code-cell}
msprime <- reticulate::import("msprime")
ts <- msprime$sim_ancestry(80, sequence_length=1e4, recombination_rate=1e-4, random_seed=42)
ts  # See "Jupyter notebook tips" below for nicer rendering
```

## Attributes and methods

`reticulate` lets us access a Python object's attributes with the `$` operator.
For example, we can access the number of samples in the tree sequence and
store it into an R object:

```{code-cell}
n <- ts$num_samples
n
is(n)  # Show the type of the object
```

The `$` operator can also be used to call a Python object's methods, for example
the {meth}`~TreeSequence.simplify` method for a tree sequence.
Method arguments are given as native R objects and then passed to `reticulate` Python
(but note that object IDs still use `tskit`'s 0-based indexing system).

```{code-cell}
reduced_ts <- ts$simplify(0:7)  # only keep samples with IDs 0, 1, 2, 3, 4, 5, 6, 7
reduced_ts <- reduced_ts$delete_intervals(list(c(6000, 10000)))  # delete data after 6 kb
reduced_ts <- reduced_ts$trim()  # remove the deleted region
paste(
    "Reduced from", ts$num_trees, "trees over", ts$sequence_length/1e3, "kb to",
    reduced_ts$num_trees, "trees over", reduced_ts$sequence_length/1e3, "kb.")
```

## IDs and indexing

If you pass a bare number to one of these methods, R treats it as a floating
point value (the default `numeric` type in R). This matters when calling
`tskit` methods that require integers (such as object IDs). For example, the
following raises a `TypeError` because `0` is passed as a float:

```{code-cell}
try(ts$node(0))  # Will raise a TypeError
is(0)
```

To pass `0` as an integer, either coerce it with `as.integer` or append `L`:

```{code-cell}
is(as.integer(0))
ts$node(as.integer(0))
# or
is(0L)
ts$node(0L)
```

You only need this coercion when calling underlying `tskit` methods that expect
integer arguments. You do not need it for array indexing itself. However, when
indexing `tskit` arrays in R, remember that `tskit` IDs start at zero, whereas
R indexes start at one:

```{code-cell}
root_id_int <- ts$first()$root
is(root_id_int)
root_id_num <- as.numeric(root_id_int)

# Using integer ID will work
paste("Root time via tskit method using integer ID:", ts$node(root_id_int)$time)
# Using numeric ID will raise TypeError
try(paste("Root time via tskit method using float ID:", ts$node(root_id_num)$time))

# When indexing into tskit arrays in R, add 1 to the ID (integer or numeric)
paste("Root time via integer array:", ts$nodes_time[root_id_int + 1])
paste("Root time via float   array:", ts$nodes_time[root_id_num + 1])
```

## Analysis

From within R, we can use `tskit`'s powerful
[Statistics](https://tskit.dev/tskit/docs/stable/stats.html) framework to
efficiently compute many summary statistics from a tree sequence. To illustrate
this, we'll first add mutations with the {func}`msprime:msprime.sim_mutations`
function and then compute genetic diversity:

```{code-cell}
ts_mut <- msprime$sim_mutations(reduced_ts, rate=1e-4, random_seed=321)
paste(ts_mut$num_mutations, "mutations, genetic diversity is", ts_mut$diversity())
```

Numerical arrays and matrices work as expected:
they are converted to equivalent R objects.
For instance, we can use the tree sequence {meth}`~TreeSequence.genotype_matrix()`
method to return the genotypes of the tree sequence as a matrix object in R.

```{code-cell}
G <- ts_mut$genotype_matrix()
G
is(G)
```

We can then use R functions directly on the genotype matrix:

```{code-cell}
allele_frequency <- rowMeans(G)
allele_frequency
```

## Jupyter notebook tips

When running R in a [Jupyter notebook](https://jupyter.org), you can define a
few functions so that `tskit` objects render nicely in notebook output:

```{code-cell}
# Define some magic functions to allow objects to be displayed in R Jupyter notebooks
repr_html.tskit.trees.TreeSequence <- function(obj, ...){obj$`_repr_html_`()}
repr_html.tskit.trees.Tree <- function(obj, ...){obj$`_repr_html_`()}
repr_svg.tskit.drawing.SVGString <- function(obj, ...){obj$`__str__`()}
```

This leads to much nicer tabular summaries:

```{code-cell}
ts_mut
```

It also allows trees and tree sequences to be plotted inline:

```{code-cell}
ts_mut$draw_svg(y_axis=TRUE, y_ticks=0:10)
```

## Interaction with R libraries

R has a number of libraries for genomic data and trees. Below we show the
phylogenetic tree representation from the popular
[ape](http://ape-package.ird.fr) package, using all trees
{meth}`exported in Nexus format<TreeSequence.write_nexus>`, or
individual trees {meth}`exported in Newick format<Tree.as_newick>`:

```{code-cell}
file <- tempfile()
ts_mut$write_nexus(file)
# Warning: ape trees are stored independently, so this uses much more memory than tskit
trees <- ape::read.nexus(file, force.multi=TRUE)  # return a set of trees

# Or simply read in a single tree
tree <- ape::read.tree(text=ts_mut$first()$as_newick())

# Now we can plot the tree in tskit style, but using the ape library
plot(tree, direction="downward", srt=90, adj=0.5)  # or equivalently use trees[[1]]
```

Note that nodes are labelled with the prefix `n`, so nodes `0`, `1`, `2`, ...
become `n0`, `n1`, `n2`, ... This helps avoidinng confusion between the
zero-based counting system used natively by `tskit` and the one-based counting
system used in `R`.

## Further information

Be sure to check out the [reticulate](https://rstudio.github.io/reticulate)
documentation, in particular on
[Calling Python from R](https://rstudio.github.io/reticulate/articles/calling_python.html),
which includes important information on how R data types are converted to their
equivalent Python types.

## Advanced usage through tskit's C API

While the `reticulate` approach is powerful and useful for most analyses, it
does add overhead from Python calls and object conversion. This overhead may be
minimal for some use cases but significant for others, such as repeatedly
calling Python from an R loop.
An alternative is to use `RcppTskit`, which provides access to the `tskit` C API
(see <https://cran.r-project.org/package=RcppTskit/>), but this is an advanced
topic beyond the scope of this tutorial.
