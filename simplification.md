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
def create_notebook_data():
    pass

# create_notebook_data()  # uncomment to recreate the tree seqs used in this notebook
```

(sec_simplification)=

# Simplification
% remove underscores in title when tutorial is complete or near-complete

**Georgia Tsambos**

:::{todo}
Create content. See https://github.com/tskit-dev/tutorials/issues/52
:::

Consider two of the most common 'subset'-type operations we might want to perform on genomic datasets:

 - Look at the data for a subset of the samples in the original dataset.
 - Look at sequence or variant information at the specific sites that vary within that subsample.

`simplify` is the tree sequence version of these operations,
but it is also more flexible than this.
Essentially, `simplify` allows you to prune away certain types of information in the tree sequence that may be irrelevant in your particular application.

```{code-cell}
import tskit
import msprime
```

## An example dataset

To demonstrate `simplify` in action, we'll simulate a scenario involving three modern-day populations, `SMALL`, `BIG` and `ADMIX`, and one ancestral population `ANC`:

```{code-cell}
demography = msprime.Demography()
demography.add_population(name="SMALL", initial_size=2000)
demography.add_population(name="BIG", initial_size=5000)
demography.add_population(name="ADMIX", initial_size=2000)
demography.add_population(name="ANC", initial_size=5000)
demography.add_admixture(
    time=100, derived="ADMIX", ancestral=["SMALL", "BIG"],
    proportions=[0.5, 0.5])
demography.add_population_split(
    time=1000, derived=["SMALL", "BIG"], ancestral="ANC")
```

Our simulated tree sequence will contain genomic information for 100 diploid individuals from each of the three contemporary populations:

```{code-cell}
ts = msprime.sim_ancestry(
    samples={"SMALL": 100, "BIG": 100, "ADMIX" : 100},
    demography=demography,
    random_seed=2432,
    sequence_length=5e7,
    recombination_rate=1e-8
)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=6151)
```
Before moving on, have a quick look at the number of elements and overall size of this tree sequence. We'll see that applying `simplify` will reduce many of these.

```{code-cell}
ts
```

## The basic syntax

At minimum, `simplify` requires a list of sample IDs that you wish to include in the new, smaller tree sequence.
Suppose we wanted to create a tree sequence holding the coalescent history of only those samples from population 'ADMIX' (which has a population label of 2):

```{code-cell}
tss = ts.simplify(ts.samples(2))
tss
```

We now have a smaller tree sequence holding just those 200 sample chromosomes from population 'ADMIX'.

Although there are fewer edges and nodes in this newer, simplified tree sequence, and the total file size is smaller,
it's not a drastic difference.
For instance, we have reduced the total number of nodes by less than half (from 48663 to 30390),
even though we are looking at just a third of our original samples.
This demonstrates the sub-linearity and efficiency of tree sequence structures.

Note that the numbers of sites and mutations are also reduced.
This is because `simplify` has removed all the mutations on the edges that were pruned away, and all of the corresponding sites (we'll see how to make `simplify` behave differently later on).
The only mutations that remain in this smaller tree sequence are those
that are inherited by *some, but not all* of the admixed samples.
That is, the mutations that produced *variation* within this subsample.

## When is `simplify` most useful?

Suppose that there are certain calculations that we only wish to perform on the subsample of admixed genomes.

:::{todo}
Think about this a bit more. Main overall 'applications' are (1) making tree sequences nicer and more regular ('simpler'), (2) pruning away irrelevant stuff so that other operations run more quickly/efficiently, and (3) forward-time sims (though this is a specific case of (2)). How much quicker/more efficient depends on the scaling of the operations -- probably a big difference for something like ibd_segments, less for tree stats stuff.
:::

## Keeping unary nodes 

## Keeping invariant sites




