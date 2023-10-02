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

(sec_tskit_more_forward_simulations)=

# _Advanced forward simulations_

% remove underscores in title when tutorial is complete or near-complete

:::{todo}
Add further details on building a forward simulator
(see issue [#14](https://github.com/tskit-dev/tutorials/issues/14))
:::

In the {ref}`previous tutorial<sec_tskit_forward_simulations>`, we developed a basic
basic forward-time Wright-Fisher (WF) simulator (refer back to that tutorial for a
detailed run through of the hidden code):

```{code-cell}
:tags: ["hide-cell"]
import tskit
import numpy as np

random_seed = 6
random = np.random.default_rng(random_seed)  # A random number generator for general use

L = 50_000  # The sequence length: 50 Kb

def add_inheritance_paths(tables, parent_genomes, child_genome, recombination_rate):
    "Add paths from parent genomes to the child genome, with crossover recombination"
    L = tables.sequence_length
    num_recombinations = random.poisson(recombination_rate * L)
    breakpoints = random.integers(0, L - 1, size=num_recombinations)
    break_pos, counts = np.unique(breakpoints, return_counts=True)
    crossovers = break_pos[counts % 2 == 1]  # no crossover if e.g. 2 breaks at same pos
    left_positions = np.insert(crossovers, 0, 0)
    right_positions = np.append(crossovers, L)

    inherit_from = random.integers(2)
    for left, right in zip(left_positions, right_positions):
        tables.edges.add_row(
            left, right, parent_genomes[inherit_from], child_genome)
        inherit_from = 1 - inherit_from  # switch to other parent genome

def make_diploid(tables, time, parent_individuals=None):
    individual_id = tables.individuals.add_row(parents=parent_individuals)
    return individual_id, (
        tables.nodes.add_row(time=time, individual=individual_id),
        tables.nodes.add_row(time=time, individual=individual_id),
    )

def new_population(tables, time, prev_pop, recombination_rate):
    pop = {}
    prev_individuals = np.array([i for i in prev_pop.keys()], dtype=np.int32)
    for _ in range(len(prev_pop)):
        mother_and_father = random.choice(prev_individuals, 2, replace=True)
        child_id, child_genomes = make_diploid(tables, time, mother_and_father)
        pop[child_id] = child_genomes  # store the genome IDs
        for child_genome, parent_individual in zip(child_genomes, mother_and_father):
            parent_genomes = prev_pop[parent_individual]
            add_inheritance_paths(tables, parent_genomes, child_genome, recombination_rate)
    return pop

def initialise_population(tables, time, size) -> dict:
    return dict(make_diploid(tables, time) for _ in range(size))

```

The main simulation function, as below, returned an unsimplified tree sequence,
which we subsequently {meth}`simplified<TreeSequence.simplify>`:

```{code-cell} ipython3

def forward_WF(num_diploids, seq_len, generations, recombination_rate=0, random_seed=7):
    global random
    random = np.random.default_rng(random_seed) 
    tables = tskit.TableCollection(seq_len)
    tables.time_units = "generations"

    pop = initialise_population(tables, generations, num_diploids)
    while generations > 0:
        generations = generations - 1
        pop = new_population(tables, generations, pop, recombination_rate)

    tables.sort()
    return tables.tree_sequence()
```

## Repeated simplification

We can perform simplification directly on the tables within the `forward_WF()` function,
using {meth}`TableCollection.simplify`. More importantly, we can carry this out at
repeated intervals. It is helpful to think of this as regular "garbage collection",
as what we're really doing is getting rid of extinct lineages while also "trimming"
extant lineages down to a minimal representation.

:::{caution}
Regular garbage collection forces us to reckon with the fact that simplification
{ref}`changes the node IDs <sec_tskit_forward_simulations_simplification_id_changes>`.
We therefore need to remap any node (and individual) IDs that are used outside of
`tskit`. In the implementation described here, those IDs are stored in the `pop`
variable.
:::

```{code-cell}
def simplify_tables(tables, samples, pop) -> dict[int, tuple[int, int]]:
    """
    Simplify the tables with respect to the given samples, returning a
    population dict in which individual and nodes have been remapped to their
    new ID numbers
    """
    tables.sort()
    node_map = tables.simplify(samples, keep_input_roots=True)
    
    nodes_individual = tables.nodes.individual
    remapped_pop = {}
    for node1, node2 in pop.values():
        node1, node2 = node_map[[node1, node2]]  # remap
        assert nodes_individual[node1] == nodes_individual[node2]  # sanity check
        remapped_pop[nodes_individual[node1]] = (node1, node2)
    return remapped_pop


def forward_WF(
    num_diploids,
    seq_len,
    generations,
    recombination_rate=0,
    simplification_interval=None,  # default to simplifying only at end
    show=None,
    random_seed=7,
):
    global random
    random = np.random.default_rng(random_seed) 
    tables = tskit.TableCollection(seq_len)
    tables.time_units = "generations"  # optional, but helpful when plotting
    if simplification_interval is None:
        simplification_interval = generations
    simplify_mod = generations % simplification_interval

    pop = initialise_population(tables, generations, num_diploids)
    while generations > 0:
        generations = generations - 1
        pop = new_population(tables, generations, pop, recombination_rate)
        if generations > 0 and generations % simplification_interval == simplify_mod:
            current_nodes = [u for nodes in pop.values() for u in nodes]
            pop = simplify_tables(tables, current_nodes, pop)
            if show:
                print("Simplified", generations, "generations before end")

    pop = simplify_tables(tables, [u for nodes in pop.values() for u in nodes], pop)
    if show:
        print("Final simplification")
    return tables.tree_sequence()

ts = forward_WF(6, L, generations=100, simplification_interval=25, show=True)
ts.draw_svg(size=(800, 200))
```

### Invariance to simplification interval
A critical concept to keep in mind is that the simulation itself is the only random component.
The simplification algorithm is deterministic given a set of (nodes, edges) satisfying
`tskit`'s sorting requirements. Therefore, the results of our new `forward_WF()` function
must be the same for all simplification intervals

:::{note}
This invariance property only holds in some cases.
We discuss this in more detail below when we add in mutation.
:::

```{code-cell}
ts = forward_WF(10, L, 500, simplification_interval=1, random_seed=42)

# Iterate over a range of odd and even simplification intervals.
print("Testing invariance to simplification interval")
test_intervals = list(range(2, 500, 33))
for i in test_intervals:
    # Make sure each new sim starts with same random seed!
    ts_test = forward_WF(10, L, 500, simplification_interval=i, show=False, random_seed=42)
    assert ts.equals(ts_test, ignore_provenance=True)
print(f"Intervals {test_intervals} passed")
```

:::{tip}
Testing your own code using loops like the one above is a very
good way to identify subtle bugs in book-keeping.
:::

### Summary

* Simplifying during a simulation changes IDs in the tree sequence tables, so we need to remap
entities that store any of these IDs between generations.
* Our code to carry out simplification gets called both during the simulation and at the end.
It's therefore worth encapsulating it into a class or function for easier code re-use and testing.

#### Technical notes

We have found that it is possible to write a simulation where the results differ
by simplification interval, but appear correct in distribution. 
By this we mean that looking at distributions of numbers of mutations, their frequencies, etc.,
match predictions from analytical theory.  However, our experience is that such simulations
contain bugs and that the summaries being used for testing are too crude to catch them.
For example, they may affect the variance in a subtle way that would require millions
of simulations to catch.  Often what is going on is that parent/offspring relationships
are not being properly recorded, resulting in lineages that either persist too long or
not long enough.  (In other words, the variance in offspring number per diploid is no
longer what it should be, meaning you've changed the effective population size.)
Thus, please make sure you get the **same** `tskit` tables out of a simulation for
any simplification interval.


## Mutations

In this section, we will add mutation to our simulation.  Mutations will occur according to the
infinitely-many sites model, which means that a new mutation cannot arise at a currently-mutated
position. $\theta = 4N\mu$ is the scaled mutation rate, and is equal to twice the expected number
of new mutations per generation.  The parameter $\mu$ is the expected number of new mutations
per gamete, per generation. Mutation positions will be uniformly distributed along the genome.

Adding mutations changes the complexity of the simulation quite a bit, because now we must
add to and simplify [site tables](sec_edge_table_definition) and
[mutation tables](sec_mutation_table_definition) instances. We might also
want to add *metadata* to the sites or mutations, recording details such as
the selection coefficient of a mutation, or the type of mutation (e.g., synonymous vs. non-synonymous).

We will write a mutation function here which we will re-use in future examples.

:::{note}
We will be treating mutations as neutral.  Doing so is odd, as one big
selling point of `tskit` is the ability to skip the tracking of neutral mutations
in forward simulations.  However, tracking neutral mutations plus metadata is the
same as tracking selected mutations and their metadata, and being able to do neat
things like put your selected mutations onto a figure of the genealogy
is one of several possible use cases.
:::

:::{todo}
The rest of this tutorial is still under construction, and needs porting from
[this workbook](https://github.com/tskit-dev/tutorials/blob/main/old-content/notebooks/wfforward.ipynb).
This will primarily deal with sites and mutations (and mutational metadata).
We could also include details on selection, if that seems sensible.

The section in that workbook on "Starting with a prior history" should be put in
the {ref}`sec_completing_forwards_simulations` tutorial.
:::