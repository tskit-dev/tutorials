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

(sec_tskit_forward_simulations)=

# Building a forward simulator

This tutorial shows the basics of creating a basic forward-time simulator that
stores the evolving genealogy as a tskit tree sequence. We will focus on the case of
diploids, in which each individual contains 2 genomes, but the concepts used here
generalize to any ploidy, if you are willing to do the book-keeping!

:::{note}
If you are simply trying to obtain a tree sequence which is
the result of a forwards-in-time simulation, this can be done by using one of the
highly capable forwards-in-time genetic simulators that already exist, such as
[SLiM](https://messerlab.org/slim/) or [fwdpy11](https://github.com/molpopgen/fwdpy11).
Documentation and tutorials for these tools exist on their respective websites. This
tutorial focusses instead on illustrating the *general principles* that lie behind such
simulators.
:::


## Definitions

Before we can make any progress, we require a few definitions.

Node
: Each node represents a genome at a point in time (often we imagine this as the "birth time"
of the genome).  It can be described by a tuple, `(id, time)`, where `id` is a unique integer,
and `time` reflects the birth time of that `id`. When generating a tree sequence, this will be
stored in a row of the [node table](sec_node_table_definition).

Individual
: A diploid individual is a group of two nodes. During simulation,a simple and efficient
way to group is to assign sequential pairs of node IDs to an individual. Although not strictly
necessary, it can also be helpful to store an individual within the tree sequence as
a row of the [individual table](sec_individual_table_definition)
(an individual can be linked to a node by storing that individual's ID
in the appropriate row of the node table). Using the existing `tskit` table
structure permits efficent storage of information such as the *geographical location*
and the *parents* of an individual (which effectively stores the entire pedigree)

Edge
: An edge reflects a transmission event between nodes.  An edge is a tuple
`(Left, Right, Parent, Child)` whose meaning is "Parent genome $P$ passed on the genomic interval
$[L, R)$ to child genome $C$". In a tree sequence this is stored in a row of the
[edge table](sec_edge_table_definition).

Time
: The time, in the discrete-time Wright-Fisher model which we will simulate,
is measured in integer generations. To match the tskit notion of time, we record
time in *generations ago*: i.e. for a simple simulation of G generations, we start
the simulation at generation $G-1$ and count down until we reach generation 0
(the current-day).

Population
: The *population* consists of $N$ diploid individuals ($2N$ nodes) at a particular
time $t$. At the start, the population will have no known ancestry, but subsequently,
each individual will be formed by choosing (at random) two parent individuals from
the population in the previous generation.


## Approach

The goal of this tutorial is to work through the book-keeping required
to generate edges, nodes, and individuals forwards in time, adding them to
the relevant `tskit` tables. The simulation will then export these tables to an
immutable [tree sequence](sec_what_is) for storage or analysis.

To aid efficiency, we will describe how to {ref}`simplify<sec_simplification>` the
tree sequence into the minimal set of nodes and edges that describe
the history of "sample" set of nodes (for example, the current-day genomes).

The model we will simulate is a discrete-time Wright-Fisher (WF) model, but the principles
can be extended to many other forward-time models. The demonstration code is written in Python,
but the logic applies to forward simulations written in any language. Note that
a more concise example of a 20-line haploid forward-time simulator is given at the
top of the {ref}`sec_completing_forwards_simulations` tutorial.

### Setup

First, we'll activate the necessary libraries and define some general parameters.
The [numpy](https://numpy.org/doc/stable/) library will be used to produce random numbers.

```{code-cell}
import tskit
import numpy as np

random_seed = 123
random = np.random.default_rng(random_seed)  # A random number generator for general use

sequence_length = 50_000  # 50 Kb
```

## Simulation without recombination

We'll start by simulating a small region of a larger genome (e.g. a "gene", or a portion of a gene).

```{code-cell}
focal_region = [20_000, 21_000]
```

The assumption here is that the region is small enough that there is no recombination.
First we will define a function `add_inheritance_paths()` that takes one of the child's
genomes in our simulation, say the maternal one, and creates `tskit` edges that
connect it to the two genomes present in the mother. With no recombination,
we can create a single edge, linked to one of the mother's genomes at random,
The function can then be applied a second time focussing on the paternal genome.

```{code-cell}
def add_inheritance_paths(tables, parent_genomes, child_genome):
    "Add inheritance paths from one of two parent genomes to the child genome"
    assert len(parent_genomes) == 2
    left, right = focal_region  # only define inheritance in this focal region
    inherit_from = random.integers(2)  # randomly chose parent 0 or 1
    tables.edges.add_row(left, right, parent_genomes[inherit_from], child_genome)
```

### Core steps

The core simulation function generates a new population by repeating the following steps:
1. Pick a pair of parent individuals at random from the population in the previous generation.
2. Create a new child individual with two genomes (an individual ID will be created
   by adding a row to the [individual table](sec_individual_table_definition)
   and two nodes IDs will be created by adding two rows to the
   [node table](sec_node_table_definition)). If desired, we can also provide
   the parent IDs when adding to the individual table, which will result in
   storing the entire genealogical pedigree.
3. Add the inheritance paths for each child genome using the `add_inheritance_paths()`
   function we defined above.

We will also define a simpler function, involving just step 2, to create the initial
population.

For convenience, the population will be stored in a Python dictionary that maps
the individual ID to the IDs of its two genomes. Note that in lower-level
languages such as C, rather than use a mapping of individual ID to two genomes,
we could instead rely on the fact that genome IDs are allocated sequentially, and
pick a random pair of genomes from the previous population directly,
without using individual IDs.

```{code-cell}
# For visualising unsimplified tree sequences, it can help to flag all nodes as samples
# (simplifying resets node flags, so setting default node_flags here is optional)
default_node_flags = tskit.NODE_IS_SAMPLE

def make_diploid(tables, time, parent_individuals=None) -> tuple[int, tuple[int, int]]:
    """
    Make an individual and its diploid genomes by adding to tables, returning the IDs.
    Specifying parent_individuals is optional but restults in the pedigree being stored.
    """
    individual_id = tables.individuals.add_row(parents=parent_individuals)
    return individual_id, (
        tables.nodes.add_row(default_node_flags, time, individual=individual_id),
        tables.nodes.add_row(default_node_flags, time, individual=individual_id),
    )

def new_population(tables, time, prev_population) -> dict[int, tuple[int, int]]:
    population = {}  # fill with individual_ID: (maternal_genome_ID, paternal_genome_ID)

    # Cache the list of individual IDs in the previous population, for efficiency
    prev_individuals = np.array([i for i in prev_population.keys()], dtype=np.int32)

    for _ in range(len(prev_population)):
        # 1. Pick two individual parent IDs at random, `replace=True` allows selfing
        mother_and_father = random.choice(prev_individuals, 2, replace=True)

        # 2. Get 1 new individual ID + 2 new node IDs
        child_id, child_genomes = make_diploid(tables, time, mother_and_father)
        population[child_id] = child_genomes  # store the genome IDs
        
        # 3. Add inheritance paths to both child genomes
        for child_genome, parent_individual in zip(child_genomes, mother_and_father):
            parent_genomes = prev_population[parent_individual]
            add_inheritance_paths(tables, parent_genomes, child_genome)
    return population

def initialise_population(tables, time, size) -> dict[int, tuple[int, int]]:
    # Just return a dictionary by repeating step 2 above
    return dict(make_diploid(tables, time) for _ in range(size))
```

:::{note}
For simplicity, the code above assumes any parent can be a mother or a father
(i.e. this is a hermaphrodite species). It also allows the same parent to be
chosen as a mother and as a father (i.e. "selfing" is allowed),
which gives simpler theoretical results. This is easy to change if required.
:::

### Simulation loop

The forward-in-time simulator simply repeats the `new_population()` routine,
replacing the old population with the new one. For efficiency reasons,
`tskit` has strict requirements for the order of edges in the edge table,
so we need to {meth}`~TableCollection.sort` the tables before we output the
final tree sequence.

```{code-cell}
def simple_diploid_sim(diploid_population_size, generations) -> tskit.TreeSequence:
    tables = tskit.TableCollection(sequence_length)
    tables.time_units = "generations"  # optional, but helpful when plotting

    population = initialise_population(tables, generations, diploid_population_size)
    while generations > 0:
        generations = generations - 1
        population = new_population(tables, generations, population)

    tables.sort()
    return tables.tree_sequence()


### Test it for two generations
ts = simple_diploid_sim(diploid_population_size=6, generations=1)
ts.draw_svg(y_axis=True, size=(800, 180))
```

It looks like it is working correctly: each row in this plot corresponds to a generation,
and all 12 genomes in the current generation
(time=0) trace back to a genome in the initial generation (time=1). However, not
all individuals in the initial generation have passed on genetic material at this
genomic position (they appear as isolated nodes at the top of the plot).

Now let's simulate for a longer time period, and set a few helpful plotting parameters.
By convention, we plot the most recent generation at the bottom of the plot
(i.e. perversely, each tree has "leaves" towards the bottom, and "roots" at the top)

```{code-cell}
ts = simple_diploid_sim(6, generations=15)

graphics_params = {
    "y_axis": True,
    "y_label": f"Time ({ts.time_units} ago)",
    "y_ticks": {i: 'Current' if i==0 else str(i) for i in range(16)},
}
ts.draw_svg(size=(1200, 350), **graphics_params)
```

This is starting to look like a real genealogy! However, we can also see that there is a
significant amount of redundancy in this representation: there are many nodes and edges
that do not contribute to the genealogy of the extant genomes. The process of
simplification can remove this.


(sec_tskit_forward_simulations_simplification)=

## Simplification

The key to efficent forward-time genealogical simulation is the process of
[simplification](sec_simplification), which can reduce much of the complexity
shown in the tree above. Typically, we want to remove all the lineages that
do not contribute to the current day genomes. We do this via the
{meth}`~tskit.TreeSequence.simplify` method, specifying that only the nodes
in the current generation are "samples".

```{code-cell}
current_day_genomes = ts.samples(time=0)
simplified_ts = ts.simplify(current_day_genomes, keep_unary=True, filter_nodes=False)
simplified_ts.draw_svg(size=(600, 300), **graphics_params)
```

(sec_tskit_forward_simulations_simplification_id_changes)=

### ID changes

We just simplified with `filter_nodes=False`, meaning that the tree sequence
retained all nodes even after simplification. However, many nodes are not longer
part of the genealogy; removing them means we can store fewer nodes
(although it will change the node IDs).

```{code-cell}
simplified_ts = ts.simplify(current_day_genomes, keep_unary=True)
simplified_ts.draw_svg(size=(600, 300), **graphics_params)
```

Note that the list of nodes passed to `simplify` (i.e. the current-day genomes)
have become the first nodes in the table, numbered from 0..11, and the
remaining nodes have been renumbered from youngest to oldest.

### Extra node removal

The `keep_unary=True` parameter meant that we kept intermediate ("unary") nodes,
even those that do not not represent branch-points in the tree. Often these are
also unneeded, and by default the `simplify` methods removes those too.
This will mean that the node IDs of older nodes will change again, and while we
we still retain information about genetic relationships between the *sample nodes*, we
are likely to lose most information about the parent-child relationships between
*individuals* (i.e. the "pedigree").

```{code-cell}
simplified_ts = ts.simplify(current_day_genomes)
simplified_ts.draw_svg(size=(400, 300), **graphics_params)
```

This is now looking much more like a "normal" genetic genealogy
(a "gene tree"), in which all the sample genomes trace back to a single
common ancestor.

## Multiple roots

If we run the simulation for fewer generations, it is not guaranteed
that the samples in the current day will share a common ancestor within
the timeframe of our simulation.

```{code-cell}
graphics_params["y_label"] = "Time"
ts = simple_diploid_sim(6, generations=5)
ts.draw_svg(size=(700, 200), **graphics_params)
```

This is clearer in the simplified version below. The genealogy, even at a single
point on the genome doesn't look quite like a normal "tree": instead it it contains
several unlinked topologies. In `tskit` terminology this is a single tree with
[multiple roots](sec_data_model_tree_roots)

```{code-cell}
simplified_ts = ts.simplify(ts.samples(time=0))
simplified_ts.draw_svg(size=(700, 200), **graphics_params)
```

When a forward-simulated tree has multiple roots, it can be useful to retain relevant
lineages all the way back to the start of the simulation, which allows prior history
to be grafted onto the tree sequence
(see {ref}`sec_completing_forwards_simulations_input_roots`).
This can be done using the `keep_input_roots` option:

```{code-cell}
simplified_ts = ts.simplify(ts.samples(time=0), keep_input_roots=True)
simplified_ts.draw_svg(size=(700, 200), **graphics_params)
```

## Simulation with recombination

It is relatively easy to modify this simulation code to allow recombination.
All we need to do is to redefine the `add_inheritance_paths()` function,
so that the child inherits a mosaic of the two genomes present in each parent.

The exact details of the mosaic will depend on the model of recombination you
wish to implement. For instance, a simple model such as that in the
{ref}`sec_completing_forwards_simulations` tutorial might assume exactly one
crossover per chromosome. A complex model might allow not just multiple
crossovers with e.g. recombination "hotspots", but also non-crossover events
such as {ref}`msprime:sec_ancestry_gene_conversion`.

Below is a redefined `add_inheritance_paths()` function of intermediate complexity,
which models recombination as a uniform Poisson process along a genomic interval.
We generate a set of "breakpoints" along the genome, then allocates edges from the
child genome to one or the other genome from the parent, with left and right positions
reflecting the breakpoint positions. Note the the recombination_rate is likely to be
low: in most species there are relatively small number of breakpoints per chromosome,
e.g. in humans around 1 or 2.

```{code-cell}
recombination_rate = 5e-7

def add_inheritance_paths(tables, parent_genomes, child_genome):
    "Add paths from parent genomes to the child genome, with crossover recombination"
    L = tables.sequence_length
    num_recombinations = random.poisson(recombination_rate * L)
    breakpoints = random.integers(0, L - 1, size=num_recombinations)
    break_pos, counts = np.unique(breakpoints, return_counts=True)
    crossovers = break_pos[counts % 2 == 1]  # no crossover if e.g. 2 breaks at same pos
    left_positions = np.insert(crossovers, 0, 0)
    right_positions = np.append(crossovers, L)

    inherit_from = random.integers(2)  # starting parental genome
    # iterate over pairs of ([0, b1], [b1, b2], [b2, b3], ... [bN, L])
    for left, right in zip(left_positions, right_positions):
        tables.edges.add_row(
            left, right, parent_genomes[inherit_from], child_genome)
        inherit_from = 1 - inherit_from  # switch to other parent genome


# Simulate a few generations, for testing
ts = simple_diploid_sim(6, generations=5)  # Now includes recombination
ts  # Show the tree sequence
```

You can see that recombination has led to more than one tree.
In fact, there are 2 "local" trees along the genome.
Here's how the full (unsimplified) genealogy looks:

```{code-cell}
graphics_params["y_label"] = f"Time ({ts.time_units} ago)"
ts.draw_svg(size=(1000, 300), **graphics_params)
```

This is rather confusing to visualise, and will get even worse if we simulate more generations.
However, {ref}`sec_tskit_forward_simulations_simplification` can be equally applied to a
recombinant genealogy, and will reduce the genealogy to something more managable for
analysis and for visualization. Below, for example we run the same simulation code
for many more generations (using a log scale for plotting), and simplify the resulting
tree sequence to the current-day samples.

```{code-cell}
ts = simple_diploid_sim(6, generations=50)
simplified_ts = ts.simplify(ts.samples(time=0), keep_input_roots=True)
graphics_params["y_ticks"] = {0: "current", 1: 1, 2: 2, 5: 5, 10: 10, 20: 20, 50: 50}
simplified_ts.draw_svg(size=(1000, 250), time_scale="log_time", **graphics_params)
```

### Larger simulations

So far we have only simulated a relatively small population size (12 genomes / 6 diploids)
for a short time. We can easily simulate for longer times, and increase the population size.
In this case, simplification
can have a dramatic effect on the disk storage and memory required:

```{code-cell}
from datetime import datetime

population_size = 100
gens = 1000  # ten times the population size

start = datetime.now()
large_ts = simple_diploid_sim(population_size, generations=500)
print(
    f"Simulated {population_size} individuals ({population_size * 2} genomes)",
    f"for {gens} generations in {(datetime.now() - start).seconds:.1f} seconds",
)
```

```{code-cell}
print(f"Full tree sequence including dead lineages: {large_ts.nbytes/1024/1024:.2f} MB")
current_day_genomes = large_ts.samples(time=0)
simplified_ts = large_ts.simplify(current_day_genomes, keep_input_roots=True)
print(
    f"Tree sequence of current-day individuals: {simplified_ts.nbytes/1024/1024:.2f} MB,",
    f"{simplified_ts.num_trees} trees."
)
print(
    "Simplification has reduced the size by a factor of",
    f"{large_ts.nbytes / simplified_ts.nbytes:.2f}"
)
```

The most obvious improvement when simulating genealogies in forward time
is therefore to carry out regular simplification steps *during* the simulation,
rather than just at the end. This is
described in the next tutorial: {ref}`sec_tskit_more_forward_simulations`.

(sec_forward_simulations_ensuring_coalescence)=

## Ensuring coalescence

Even though the simulation above was run for hundreds of generations, there are still
trees with multiple roots in some regions of the genome. 1000 generations is not therefore
long enough to capture the ancestry back to a single common ancestor
(i.e. to ensure "full coalescence" of all local trees):

```{code-cell}
from matplotlib import pyplot as plt
plt.figure(figsize=(7, 2))
plt.stairs(
    [tree.num_roots for tree in simplified_ts.trees()],
    simplified_ts.breakpoints(as_array=True),
    baseline=None,
)
plt.xlabel("Genome position (bp)")
plt.ylabel("Number of roots")
plt.show()
```

Most of our trees have one root, but a handful have two. 
What is happening is that even though we simulated for $10N$ generations,
and the theoretical expectation for the time to the most recent common ancestor
is $E[tMRCA]=4N$, the variance is pretty big, such that some marginal
trees are not completely coalesced.  Let's take a look at the distribution of
tMRCAs using `msprime` to carry out an equivalent backward-time simulation
(which is much faster, so we can replicate it, say, 1000 times).

```{code-cell}
import msprime

tmrca=[]
sequence_length = 50_000
recombination_rate=5e-7
for sim in msprime.sim_ancestry(
    100,
    sequence_length=sequence_length,
    recombination_rate=recombination_rate,
    random_seed=42,
    num_replicates=1000,
):
    for tree in sim.trees():
        tmrca.append(tree.time(tree.root))

plt.figure(figsize=(7, 3))
plt.hist(tmrca, bins=30, density=True)
plt.xlabel("tMRCA (units of N generations)")
plt.ylabel("Proportion of trees")
plt.show()
```

So, some fraction of local trees have very large tMRCAs! And of course
the larger the population, the longer the time needed to ensure full coalescence.
Hence for large forward simulation models, the number of generations required 
to ensure full coalescence can be prohibitive.

Uncoalesced local trees matter, because they result in some pairs of
samples which are not connected, and so the genetic distance betweem
them is undefined (another way to imagine this is that we cannot spot
mutations that occur above the local roots of the tree). Hence the
simulation is not capturing the full genetic diversity within the sample.

A powerful way to get around this problem is *recapitation*,
in which an alternative technique, such as backward-in-time coalescent simulation,
is used to to fill in the "head" of the tree sequence. In other words,
we can use a fast backward-time simulator such as `msprime` to simulate the
prior genealogy of the oldest nodes in the simplified tree sequence.
Details are described in the {ref}`sec_completing_forwards_simulations`
tutorial.

## More complex forward-simulations

The next tutorial shows the principles
behind more complex simulations, e.g. regular simplification during simulation,
adding mutations, and adding metadata. It also details several
extra tips and tricks we have learned when building forward simulators.