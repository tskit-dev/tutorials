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

def computing_statistics():
    ts = msprime.simulate(
        10**4, Ne=10**4, recombination_rate=1e-8, mutation_rate=1e-8, length=10**7, random_seed=42)
    ts.dump("data/computing_statistics.trees")
    
def afs():
    ts = msprime.simulate(6, mutation_rate=1, random_seed=47)
    # remove the mutation times so the plot is nicer
    tables = ts.dump_tables()
    tables.mutations.time = np.full_like(tables.mutations.time, tskit.UNKNOWN_TIME)
    ts = tables.tree_sequence()
    ts.dump("data/afs.trees")

def create_notebook_data():
    computing_statistics()
    afs()

# create_notebook_data()  # uncomment to recreate the tree seqs used in this notebook
```

(sec_analysing_tree_sequences)=

# _Analysing tree sequences_
% remove underscores in title when tutorial is complete or near-complete

:::{note}
This tutorial is a work in progress.
:::

(sec_tutorial_stats)=

## Computing statistics

Tskit provides an extensive and flexible interface for computing population
genetic statistics, which is documented in detail in the
{ref}`general statistics<tskit:sec_stats>` section of the offical documentation.
This tutorial aims to give a quick overview of how the statistics APIs work
and how to use them effectively. After this overview, we will dive into
allele frequency spectrum analysis.

First, let's load a tree sequence that represents 10 thousand samples
over a 10 megabases chromosome from a simulation with human-like parameters:

```{code-cell} ipython3
import tskit
ts = tskit.load("data/computing_statistics.trees")
print(ts)
```

This tree sequence has ~36.6 thousand trees & ~39 thousand segregating sites.
We'd now like to analyse this dataset using different statistics.

### One-way statistics

We refer to statistics that are defined with respect to a single set of
samples as "one-way". An example of such a statistic is diversity, which
is computed using the {meth}`TreeSequence.diversity` method:

```{code-cell} ipython3
d = ts.diversity()
print(f"Average diversity per unit sequence length = {d:e}")
```

The above result tells us the average diversity across the samples and
their whole sequence (that is their genomes), returning a single number.
Note the emphasis on "average" in the above print statement.
See the {ref}`Span normalisation<sec_span_normalisation>` sub-section
below for more details on this point.

#### Windows

We often want to compute statistics in {ref}`windows<tskit:sec_stats_windows>`
along the genome. We can use the ``windows`` argument to do this:

```{code-cell} ipython3
import numpy as np
windows = np.linspace(0, ts.sequence_length, num=5)
print(windows)
d = ts.diversity(windows=windows)
print(d)
```

The ``windows`` argument takes a numpy array with breakpoints along the genome.
Above, we created four equally spaced windows of size 2.5 megabases,
since the total sequence length (``ts.sequence_length``) is 10 megabases.
Note that the array contains `k + 1` elements to define `k` windows.
Because we have asked for values in windows, tskit now returns a numpy array
rather than a single value.
See {ref}`tskit:sec_stats_windows` for a full description on the ``windows`` argument,
including shortcuts on how to get statistics per local tree or per site.
Note though that there can be many local trees and sites in a tree sequence,
leading to a lot of output.
See {ref}`tskit:sec_stats_output_dimensions` for a full description of how
the output dimensions of statistics are determined by the ``windows`` argument.

#### Sample sets

Suppose we wanted to compute average diversity within a specific subset of samples,
instead of all samples.
We can do this using the {ref}`sample_sets<tskit:sec_stats_sample_sets>` argument:

```{code-cell} ipython3
A = ts.samples()[:100]
d = ts.diversity(sample_sets=A)
print(d)
```

Here, we've computed the average diversity within the first hundred samples across
the whole sequence. As we've not specified any windows, this is again a single value.

We can also compute average diversity in *multiple* sample sets at the same time
by providing a list of sample sets as an argument:

```{code-cell} ipython3
A = ts.samples()[:100]
B = ts.samples()[100:200]
C = ts.samples()[200:300]
d = ts.diversity(sample_sets=[A, B, C])
print(d)
```

Because we've computed a statistic for multiple sets, tskit returns a numpy array.
We have asked for diversity within three different sample sets,
and tskit therefore returns an array with three values.
As with ``windows``, the ``sample_sets`` argument determines the output dimensions.
See {ref}`tskit:sec_stats_output_dimensions` for a detailed description of the rules.

#### Sample sets and windows

We can also compute a statistic for multiple sets in multiple windows:

```{code-cell} ipython3
d = ts.diversity(sample_sets=[A, B, C], windows=windows)
print("shape = ", d.shape, "\n", d)
```

We have computed average diversity in three sample sets across four genomic windows,
and our output is therefore a 2D numpy array with four rows and three columns:
each row contains the diversity values within sample set A, B, and C
for a particular window.

(sec_span_normalisation)=

#### Span normalisation

Above we referred to the *average* diversity per unit sequence length,
either for the whole sequence or within windows.
To facilitate comparison, tskit by default normalises statistics
by the sequence length, returning average of a statistic.
When this is not what you want,
you can use the {ref}`span_normalise<tskit:sec_stats_span_normalise>` argument:

```{code-cell} ipython3
davg = ts.diversity()
dtot = ts.diversity(span_normalise=False)
print(f"Average diversity per unit sequence length = {davg:e}")
print(f"Total diversity across sequence = {dtot}")
print(davg * ts.sequence_length)
```

As you see above, the average diversity per unit sequence length is equal to
the total diversity divided (normalised) by the sequence length.

Normalisation is useful when comparing statistics across different windows
and/or sample sets. As an example, let's look at the span of
the first four local trees and compute the total diversity for each tree:

```{code-cell} ipython3
b = ts.breakpoints(as_array=True)[0:4]
print(b)
dtot = ts.diversity(windows="trees", span_normalise=False)[0:4]
print(f"Total diversity per local tree = {dtot}")
```

As you can see the span of local tree differs.
The first tree spans from 0 to 131.2...,
the second tree spans from 131.2... to 147.8...,
etc.
Consequently, the total diversity for each tree is also different.
The first two trees have total diversity of 0, indicating they contain no mutations.
The third tree has total diversity of 0.225... and
the fourth tree has total diversity of 0.001...
These last two values differ because of diffent tree spans and frequency of mutated
alleles in these trees, challenging comparison of the values.
To facilitate comparisons, tskit by default span-normalises statistics:

```{code-cell} ipython3
davg = ts.diversity(windows="trees")[0:4]
print(f"Average diversity per local tree = {davg}")
```

The above values indicate that the average diversity per unit sequence length
is larger in the third tree than in the fourth tree, indicating that the third tree
has more intermediate frequencies of alleles than the fourth tree.

#### Mode

TODO: The simulated tree sequence is without time units,
hence branch statistics have no units. Should we skip
the example below then and just keep the 1st paragraph below
(but surely someone will try branch mode and ask what is the meaning of
ts.diversity(mode="branch") being 40534.4...)?

Above we have computed statistics based on site mutations in the tree sequence.
In tskit, these statistics are hence called "site-based". This "mode" of statistics
is computed by default, but we can also compute "branch-based" and "node-based"
statistics using the {ref}`mode<tskit:sec_stats_mode>` argument:

```{code-cell} ipython3
ds = ts.diversity()
print(ds)
db = ts.diversity(mode="branch")
print(db)
```

Site-based statistics use mutation differences between genomes,
while branch-based statistics use branch lengths between genomes.
Above, 0.0004... is average diversity per unit sequence length
in terms of mutation differences between samples. And, 40534.4...
is average diversity per unit sequence length in terms of
branch lengths between samples.

Node-based statistics behave differently. They are computed per
each node in a tree sequence and TODO: REPORT WHAT?
Since there can be many nodes in a tree sequence, the output
can be large.

```{code-cell} ipython3
dn = ts.diversity(mode="node")
print(dn)
```

### Multi-way statistics

Many population genetic statistics compare sample sets to each other,
leading to "multi-way" statistics.
For example, the {meth}`TreeSequence.divergence` method computes
the divergence between two sets of samples:

```{code-cell} ipython3
A = ts.samples()[:100]
B = ts.samples()[100:200]
d = ts.divergence(sample_sets=[A, B])
print(d)
```

The divergence between two sets of samples A and B is a single number,
hence we obtain a single value as the result.

We can also compute "multi-way" statistics in windows along the genome, as before:

```{code-cell} ipython3
windows = np.linspace(0, ts.sequence_length, num=5)
d = ts.divergence(sample_sets=[A, B], windows=windows)
print(d)
```

Again, as we have defined four genomic windows along the sequence, the result is
numpy array with four values, representing the average divergence between sample
set A and B in each window.

A powerful feature of tskit's statistics API is that we can compute
multi-way statistics between multiple sets of samples simultaneously
using the {ref}`indexes<tskit:sec_stats_sample_sets_multi_way>` argument:

```{code-cell} ipython3
d = ts.divergence(sample_sets=[A, B, C], indexes=[(0, 1), (0, 2)])
print(d)
```

Here, we've specified three sample sets A, B and C, and we've computed the
divergences between A and B, and between A and C giving a numpy array of length 2.
The ``indexes`` argument is used to specify which pairs of sets we are interested in.

As before, we can combine computing multi-way statistics with multiple windows
to return a 2D numpy array:

```{code-cell} ipython3
d = ts.divergence(sample_sets=[A, B, C], indexes=[(0, 1), (0, 2)], windows=windows)
print(d)
```

Each row again corresponds to a window, which contains the average divergence
between the chosen sample sets.

If the ``indexes`` parameter is 1D array, we interpret this as specifying
a single statistic and remove the empty outer dimension:

```{code-cell} ipython3
d = ts.divergence([A, B, C], indexes=(0, 1))
print(d)
```

It's important to note that we don't **have** to remove empty dimensions: tskit
will only do this if you explicitly ask it to. Here, for example, we can keep the
output as an array with one value if we wish:

```
d = ts.divergence([A, B, C], indexes=[(0, 1)])
print(d)
```

See {ref}`tskit:sec_stats_output_dimensions` for a full description of how
the output dimensions of statistics are determined by the ``indexes`` argument.

(sec_tutorial_afs)=

## Allele frequency spectra

The allele frequency spectrum (AFS) is a fundamental tool in population genetics, and
tskit provides a flexible and powerful approach to computing such spectra.
Suppose we have the following simple tree sequence with just one local tree:

```{code-cell} ipython3
from IPython.display import display
ts = tskit.load("data/afs.trees")
tree = ts.first()
display(tree.draw_svg())
ts.tables.mutations
```

Computing the allele frequency spectrum is then achived using the
{meth}`~TreeSequence.allele_frequency_spectrum` method:

```{code-cell} ipython3
afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
print(afs)
```

This tells us that we have two singletons, six doubletons, one 3-ton, and
one 4-ton. We can see these in the above tree display:
the singletons are nodes 5 (for mutation 5) and 4 (for mutation 0),
the doubletons are nodes 0 and 1 (for mutations 2, 3, 4, 7, 8, and 9),
and so onwards.
To see where these counts come from, we can look at the above mutations table.
For example, the mutation 0 appeared at site 0 above sample node 4,
producing one of the singeltons.
Another example, the mutation 1 appeared at site 1 above node 8,
which is ancestor to sample nodes 2, 3, and 5, producing the 3-ton.

Note that the first element of the returned AFS array does *not* correspond
to the singletons ({ref}`see below for why<sec_tutorial_afs_zeroth_entry>`).
Because we know the ancestral and derived states at these sites,
we have set ``polarised`` argument to True. We can get the "folded" AFS by
setting ``polarised`` argument to False.
Because we want simple counts here and not averaged values, we set
``span_normalise=False``. As noted {ref}`above<sec_span_normalisation>`,
by default, statistics are divided by the sequence length,
so they are comparable between samples and windows, which is not what we want here.

We can also obtain AFS in windows along the genome:

```{code-cell} ipython3
afs = ts.allele_frequency_spectrum(windows=[0, 0.5, 1], polarised=True, span_normalise=False)
print(afs)
```

This time, we've asked for the number of sites at each frequency in two
equally long windows. Now we can see that in the first half of the sequence we
have one singleton, one doubleton, and one tripleton. The second half of the
sequence has one singleton, five doubletons, and one 4-ton.

### Joint spectra

We can also compute allele frequencies within multiple sets of samples,
giving the *joint allele frequency spectra*.

```{code-cell} ipython3
node_colours = {0: "blue", 2: "blue", 3: "blue", 1: "green", 4: "green", 5: "green"}
styles = [f".n{k} > .sym {{fill: {v}}}" for k, v in node_colours.items()]
tree.draw_svg(style = "".join(styles))
```

Here we've marked the samples as either blue or green (we can imagine
these belonging to different populations, for example). We can then compute
the joint AFS based on these two sets:

```{code-cell} ipython3
afs = ts.allele_frequency_spectrum(sample_sets=[[0, 2, 3], [1, 4, 5]], polarised=True)
print(afs)
```

This gave us a 2D numpy array, where each dimension corresponds to one sample set.
We see that there are:
two singletons in the second sample set and not present in the other sample set ($afs[0,1]$),
six singletons in both sets ($afs[1,1]$),
one doubleton in the first sample set that is singleton in the second sample set ($afs[1,2]$),
and one doubleton in both sets ($afs[2,2]$).

### Branch length spectra

Instead of summarising the number of sites that occur at different frequencies,
we can also compute the total branch lengths subtending a given number of samples
by setting ``mode="branch"``:

```{code-cell} ipython3
afs = ts.allele_frequency_spectrum(mode="branch", polarised=True, span_normalise=False)
print(afs)
```

Thus, the total branch length over one sample is 4.86, over two is 5.39, and so on.

(sec_tutorial_afs_zeroth_entry)=

### Zeroth and final entries in the AFS

The zeroth element of the AFS is significant when we are working with
sample sets that are a subset of all samples in the tree sequence.
For example, in the following we compute the AFS within the sample set
$[0, 1, 2]$:

```{code-cell} ipython3
afs = ts.allele_frequency_spectrum(sample_sets=[[0, 1, 2]], mode="branch", polarised=True)
print(afs)
```

Thus, the total branch length over samples 0, 1, and 2 is 5.3, and over pairs
from this set is 5.25. What does the zeroth value of 4.33 mean?
This is the total branch length over all samples that are **not** in this sample set.
By including this value, we maintain the property that for each tree,
the sum of the AFS for any sample set is always equal to the total branch length.
For example, here we compute:

```{code-cell} ipython3
print("sum afs          = ", np.sum(afs))
print("total branch len = ", tree.total_branch_length)
```

The final entry of the AFS is similar: it counts alleles (for mode="site") or
branches (for mode="branch") that are ancestral to all of the given sample set,
but are still polymorphic in the entire set of samples of the tree sequence.
Note, however, that alleles fixed among all the samples, e.g., ones above
the root of the tree, will not be included.
