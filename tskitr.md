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

To interface with `tskit` in R, we can use the [reticulate](https://rstudio.github.io/reticulate) R package, which lets you call `tskit`'s Python API within an R session. In this tutorial, we'll go through a couple of examples to show you how to get started. You'll need to install `reticulate` in your R session via `install.packages("reticulate")` and at a minimum state the required Python packages via `reticulate::py_require(c("tskit", "msprime"))`. This setting will use `reticulate`'s isolated Python virtual environment, but you can also use an alternative Python environment as described at https://rstudio.github.io/reticulate.

::::{margin}
:::{note}
The [slendr](https://www.slendr.net) R package uses `reticulate` extensively to work with tree sequences in R and also provides R wrappers for some `tskit` functions. 
:::
::::

We'll begin by simulating a small tree sequence using `msprime`.

```{code-cell}
msprime <- reticulate::import("msprime")
ts <- msprime$sim_ancestry(80, sequence_length=1e4, recombination_rate=1e-4, random_seed=42)
ts  # See "Jupyter notebook tips", below for how to render this nicely
```

## Attributes and methods

`reticulate` allows us to access a Python object's attributes via
the `$` operator. For example, we can access (and assign to a native R object)
the number of samples in the tree sequence:

```{code-cell}
n <- ts$num_samples
n
is(n)  # Show the type of the object
```

The `$` operator can also be used to call methods, for example, the
{meth}`~TreeSequence.simplify` method associated with the tree sequence.
The method parameters are given as native R objects
(but note that object IDs still use tskit's 0-based indexing system).

```{code-cell}
reduced_ts <- ts$simplify(0:7)  # only keep samples with ids 0, 1, 2, 3, 4, 5, 6, 7
reduced_ts <- reduced_ts$delete_intervals(list(c(6000, 10000)))  # delete data after 6kb
reduced_ts <- reduced_ts$trim()  # remove the deleted region
paste(
    "Reduced from", ts$num_trees, "trees over", ts$sequence_length/1e3, "kb to",
    reduced_ts$num_trees, "trees over", reduced_ts$sequence_length/1e3, "kb.")
```

## IDs and indexes

Note that if a bare digit is provided to one of these methods, it will be treated as a
floating point number (numeric in R). This is useful to know when calling `tskit` methods that
require integers (such as object IDs). For example, the following will not work:

```{code-cell}
:tags: [raises-exception, remove-output]
ts$node(0)  # Will raise a TypeError
is(0)
```

In this case, to force the `0` to be passed as an integer, you can either coerce it
using `as.integer` or simply append the letter `L`:

```{code-cell}
is(as.integer(0))
ts$node(as.integer(0))
# or
is(0L)
ts$node(0L)
```

Coercing in this way is only necessary when passing parameters to those underlying
`tskit` methods that expect integers. It is not needed to index into numeric arrays.
_However_, when using arrays, very careful attention must be paid to the fact that
`tskit` IDs start at zero, whereas R indexes start at one:

```{code-cell}
root_id_int <- ts$first()$root
is(root_id_int)
root_id_num <- as.numeric(root_id_int)

# Using integer ID will work
paste("Root time via tskit method:", ts$node(root_id_int)$time)
# Using numeric ID will not work
# paste("Root time via tskit method:", ts$node(root_id_num)$time)

# When indexing into tskit arrays in R, add 1 to the ID (integer or numeric)
paste("Root time via array access:", ts$nodes_time[root_id_int + 1])
paste("Root time via array access:", ts$nodes_time[root_id_num + 1])
```

## Analysis

From within R we can use `tskit`'s powerful
[Statistics](https://tskit.dev/tskit/docs/stable/stats.html) framework to efficiently
compute many different summary statistics from a tree sequence. To illustrate this,
we'll first add some mutations to our tree sequence with the
{func}`msprime:msprime.sim_mutations` function, and then compute the genetic diversity
for each of the tree sequence's sample nodes:

```{code-cell}
ts_mut <- msprime$sim_mutations(reduced_ts, rate=1e-4, random_seed=321)
paste(ts_mut$num_mutations, "mutations, genetic diversity is", ts_mut$diversity())
```

Numerical arrays and matrices work as expected. For instance, we can use the tree
sequence {meth}`~TreeSequence.genotype_matrix()` method to return the genotypes of
the tree sequence as a matrix object in R.

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

When running R within a [Jupyter notebook](https://jupyter.org), a few magic functions
can be defined that allow tskit objects to be rendered within the notebook:

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

R has a number of libraries to deal with genomic data and trees. Below we demo the
phylogenetic tree representation defined in the the popular
[ape](http://ape-package.ird.fr) package, taking all the trees
{meth}`exported in Nexus format<TreeSequence.write_nexus>`, or
individual trees {meth}`exported in Newick format<Tree.as_newick>`:

```{code-cell}
file <- tempfile()
ts_mut$write_nexus(file)
# Warning - ape trees are stored independently, so this will use much more memory than tskit
trees <- ape::read.nexus(file, force.multi=TRUE)  # return a set of trees

# Or simply read in a single tree
tree <- ape::read.tree(text=ts_mut$first()$as_newick())

# Now we can plot the tree in tskit style, but using the ape library
plot(tree, direction="downward", srt=90, adj=0.5)  # or equivalently use trees[[1]]
```

Note that nodes are labelled with the prefix `n`, so that nodes `0`, `1`, `2`, ...
become `n0`, `n1`, `n2`, ... This helps to avoid
confusion between the the zero-based counting system used natively
by `tskit`, and the one-based counting system used in `R`.

## Further information

Be sure to check out the [reticulate](https://rstudio.github.io/reticulate)
documentation, in particular on
[Calling Python from R](https://rstudio.github.io/reticulate/articles/calling_python.html),
which includes important information on how R data types are converted to their
equivalent Python types.

## Advanced usage through the tskit's C API

While the approach with `reticulate` is very powerful and will be useful for
most analyses, it does have an overhead due to calling Python and conversion of objects.
This overhead will be minimal for some use cases but could be significant for others.
An alternative approach is to use the `tskit` C API from R,
via the standard [R Extensions](https://cran.r-project.org/doc/manuals/R-exts.html) or
via [Rcpp](https://www.rcpp.org),
but these are advanced topics beyond the scope of this tutorial.
