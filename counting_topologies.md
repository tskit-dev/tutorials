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

(sec_counting_topologies)=

```{code-cell} ipython3
:tags: [remove-cell]
import msprime
import stdpopsim

def topologies_sim_speciestree():
    newick_species_tree = "((A:100.0,B:100.0):100.0,C:200.0)"
    demography = msprime.Demography.from_species_tree(newick_species_tree, initial_size=100)
    ts = msprime.sim_ancestry({0: 2, 1: 2, 2: 2}, demography=demography, random_seed=321)
    ts.dump("data/topologies_sim_speciestree.trees")

def topologies_sim_stdpopsim():
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model("OutOfAfrica_3G09")
    contig = species.get_contig("chr1", length_multiplier=0.0002, mutation_rate=model.mutation_rate)
    samples = {"YRI": 1000, "CEU": 1000, "CHB": 1000}
    engine = stdpopsim.get_engine("msprime")
    ts = engine.simulate(model, contig, samples, seed=321)
    ts.dump("data/topologies_sim_stdpopsim.trees")


def create_notebook_data():
    topologies_sim_speciestree()
    topologies_sim_stdpopsim()

# create_notebook_data()  # uncomment to recreate the tree seqs used in this notebook
```

# Counting topologies

**Yan Wong**

This tutorial is intended to be a gentle introduction to the combinatorial
treatment of tree topologies in `tskit`. For a more formal introduction,
see the {ref}`sec_combinatorics` section of the
[official `tskit` documentation](tskit:sec_introduction).

The *topology* of a single tree is the term used to describe the branching pattern,
regardless of the lengths of the branches. For example, both trees below have the
same topology, although the branch lengths differ:

```{code-cell}
import tskit
node_labels = {0: "a", 1: "b", 2: "c"}  # avoid confusion by using letters to label tips
tree = tskit.Tree.generate_comb(3)
display(tree.draw_svg(node_labels=node_labels, y_axis=True))

deep_tree = tskit.Tree.generate_comb(10).tree_sequence.simplify([0, 1, 2]).first()
display(deep_tree.draw_svg(node_labels=node_labels, y_axis=True))
```

:::{note}
The treatment of topologies in `tskit` is restricted to trees with a single defined root,
without nodes with a single child (i.e. trees must consist of nodes that are either leaves,
or internal nodes with two or more children).  For convenience in the examples
below, trees are drawn with the tips flagged as samples, although whether a node is a sample or
not does not change the topology of the tree.
:::

## Tree labellings and shapes

The topology of a tree also takes into account the labelling of tips, so that
the trees below, although they have the same *shape*, count as three
different topologies:

```{code-cell}
:tags: [hide-input]
from string import ascii_lowercase
from IPython.display import SVG

def str_none(s, prefix=None):
    if s is not None:
        if prefix is None:
            return str(s)
        else:
            return prefix + " = " + str(s)
    return None

def draw_svg_trees(trees, node_labels={}, x_lab_attr=None, width=100, height=150, space=10):
    w = width + space
    h = height + space
    trees = list(trees)
    s = f'<svg height="{h}" width="{w * len(trees)}" xmlns="http://www.w3.org/2000/svg">'
    s += f'<style>.x-axis {{transform: translateY({space}px)}}</style>'
    for i, tree in enumerate(trees):
        s += tree.draw_svg(
            size=(width, height),
            canvas_size=(w, h),
            root_svg_attributes={"x": i * w},
            node_labels=node_labels,
            x_label=str_none(getattr(tree.rank(), x_lab_attr or "", None), x_lab_attr)
        )
    s += '</svg>'
    return SVG(s)

draw_svg_trees(tskit.all_tree_labellings(tree), node_labels={u: ascii_lowercase[u] for u in tree.samples()})
```

These are, in fact, the only possible three labellings for a three-tip tree of that shape.
There is only one other possible shape for a three-tip tree, and for this shape,
all labelling orders are equivalent (in other words, there is only one
possible labelling):

```{code-cell}
:tags: [hide-input]
tskit.Tree.generate_star(3).draw_svg(node_labels={})
```

A 3-tip tree therefore has only four possible topologies.
These can be generated with the {func}`~tskit.all_trees` function.

```{code-cell}
generated_trees = tskit.all_trees(3)
print("For a three-tip tree there are", len(list(generated_trees)), "labelled topologies.")
```

Here they are, plotted out with their shapes enumerated from zero:

```{code-cell}
:tags: [hide-input]
draw_svg_trees(
    tskit.all_trees(3),
    node_labels={u: ascii_lowercase[u] for u in tree.samples()},
    x_lab_attr="shape"
)
```

### Enumerating shapes and labellings

For a tree with four tips, more topologies and shapes are possible. As before, we can generate the
topologies using {func}`~tskit.all_trees`. Alternatively, if we only want the (unlabelled) shapes,
we can use the {func}`~tskit.all_tree_shapes` function:

```{code-cell}
print("For a four-tip tree there are", len(list(tskit.all_trees(4))), "labelled topologies.")

generated_trees = tskit.all_tree_shapes(4)
print("These can be categorised into", len(list(generated_trees)), "shapes.")
```

Again, we can give each shape a number or *index*, starting from zero:

```{code-cell}
:tags: [hide-input]
draw_svg_trees(tskit.all_tree_shapes(4), x_lab_attr="shape")
```

Each of these shapes will have a separate number of possible labellings, and trees with
these labellings can be created using {func}`~tskit.all_tree_labellings`:

```{code-cell}
for shape_index, tree in enumerate(tskit.all_tree_shapes(4)):
    labellings = tskit.all_tree_labellings(tree)
    num_labellings = len(list(labellings))
    print(
        f"Tree shape {shape_index} for a four-tip tree has "
        f"{num_labellings} labelling{'' if num_labellings==1 else 's'}."
    )
```

Any tree topology for a tree of $N$ tips can therefore be described by a
shape index combined with a labelling index. This is known as the
*rank* of a tree, and it can be obtained using the
{meth}`Tree.rank` method. For instance, here is the rank of a simulated tree
of 10 tips:

```{code-cell}
:tags: [hide-input]
import msprime
num_tips = 10
simulated_ts = msprime.sim_ancestry(10, ploidy=1, random_seed=123)
simulated_tree = simulated_ts.first()
print("The topology of the simulated tree below can be described as", simulated_tree.rank())
ascii_node_labels = {u: ascii_lowercase[u] for u in simulated_tree.samples()}
simulated_tree.draw_svg(node_labels=ascii_node_labels)
```


A tree with the same topology (i.e. the same shape and labelling, but ignoring
the branch lengths) can be generated using the {meth}`Tree.unrank` method, by
specifying the number of tips and the appropriate `(shape, labelling)` tuple:

```{code-cell}
new_tree = tskit.Tree.unrank(num_tips, (1270, 21580))
new_tree.draw_svg(node_labels=ascii_node_labels)
```

Note that this method generates a single tree in a new tree sequence
whose a default sequence length is 1.0.

## Methods for large trees

The number of possible topologies for a tree with $N$ tips
grows very rapidly with $N$. For instance, with 10 tips, there are
282,137,824 possible topologies.

For this reason, the {func}`~tskit.all_trees`, {func}`~tskit.all_tree_shapes` and
{func}`~tskit.all_tree_labellings` methods do not return a list of trees
but an iterator over the trees. This means it is perfectly possible to start
iterating over (say) all tree shapes for a tree of 100 leaves, but
the iterator will not finish before the death of our galaxy.

```{code-cell}
for num_trees, tree in enumerate(tskit.all_tree_shapes(100)):
    shape = tree.rank().shape
    b2 = tree.b2_index()
    print(f"A 100-tip tree with shape index {shape} has a b2 balance index of {b2}")
    if num_trees > 5:
      break  # better not let this run too long!
```

For similar combinatorial reasons, the {meth}`Tree.rank` method can be
inefficient for large trees. To compare the topology of two trees, you are
therefore recommended to use e.g. the {meth}`Tree.kc_distance` method
rather than comparing ranks directly.

```{code-cell}
simulated_tree = simulated_ts.first(sample_lists=True)  # kc_distance requires sample lists
if simulated_ts.first(sample_lists=True).kc_distance(simulated_tree) == 0:
    print("Trees are identical")
    # To compare to the new_tree we need to fix 
    # print("The simulated and topology-constructed trees have the same topology")
```

Despite the combinatorial explosion associated with topologies of
many-tip trees, it is still possible to efficiently count
the number of *embedded topologies* in a large tree.

### Embedded topologies

An embedded topology is a a topology involving a subset of the tips of a tree.
If the tips are classified into (say) three groups, red, green, and blue,
we can efficiently count all the embedded three-tip trees which have
one tip from each group using the {meth}`Tree.count_topologies` method.

```{code-cell}
:tags: [hide-input]
big_tree = tskit.load("data/topologies_sim_speciestree.trees").first()
# Check all observed topologies have the same counts
assert list(big_tree.count_topologies()[0, 1, 2].values()) == [32, 32]
styles = [
    f".node.sample.p{p.id} > .sym " + "{" + f"fill: {colour}" + "}"
    for colour, p in zip(['red', 'green', 'blue'], big_tree.tree_sequence.populations())
]
big_tree.draw_svg(style="".join(styles), node_labels={}, time_scale="rank", x_label="big_tree")
```

In this tree, it is clear that the green and blue tips never cluster together.
The {meth}`Tree.count_topologies` method exhaustively looks at all
combinations of one red, one blue, and one green tip, and confirms that we never see
the topology grouping green and blue. However, as might be expected from
examination of the plot above, a red tip is equally likely to be a sister to a
green tip as to a blue tip:

```{code-cell}
# By default `count_topologies` chooses one tip from each population, like setting
# sample_sets=[ts.samples(p.id) for p in ts.populations() if len(ts.samples(p.id)) > 0]

topology_counter = big_tree.count_topologies()

colours = ['red', 'green', 'blue']
styles = [f".n{u}>.sym {{fill: {c} }}" for u, c in enumerate(colours)]

embedded_counts = topology_counter[0, 1, 2]
for embedded_tree in tskit.all_trees(3):
    rank = embedded_tree.rank()
    number_of_instances = embedded_counts[rank]
    label = f"{number_of_instances} instances embedded in big_tree"
    display(embedded_tree.draw_svg(style="".join(styles), node_labels={}, x_label=label))
```

## Methods over tree sequences

It can be useful to count embedded topologies over an entire tree sequence.
For instance, we might want to know the number of embedded topologies
that support Neanderthals as a sister group to europeans versus africans.
`Tskit` provides the efficient {meth}`TreeSequence.count_topologies` method to
do this [incrementally](sec_incremental), without having to re-count the topologies
independently in each tree.

```{code-cell}
:tags: [hide-input]
from myst_nb import glue
ts = tskit.load("data/topologies_sim_stdpopsim.trees")
print(f"Loaded a stdpopsim of {ts.num_trees} African+European+Chinese trees, each with {ts.num_samples} tips")
glue("seq_len", int(ts.sequence_length/1000), display=False)
```

Although the trees in this tree sequence are very large, counting the embedded topologies is
quite doable (for speed in this demo we are only simulating {glue:}`seq_len` kilobases, but
calculating the average over an entire chromosome simply takes a little longer)

```{code-cell}
from datetime import datetime
names = {"YRI": "African", "CEU": "European", "CHB": "Chinese"}
colours = {"YRI": "yellow", "CEU": "green", "CHB": "blue"}

population_map = {p.metadata["id"]: p.id for p in ts.populations()}
sample_populations = list(sorted({ts.node(u).population for u in ts.samples()}))
topology_span = {tree.rank(): 0 for tree in tskit.all_trees(len(sample_populations))}

start = datetime.now()
total = 0
for topology_counter, tree in zip(ts.count_topologies(), ts.trees()):
    embedded_topologies = topology_counter[sample_populations]
    weight = tree.span / ts.sequence_length
    for rank, count in embedded_topologies.items():
        topology_span[rank] += count * weight
        total += count
print(f"Counted {total} embedded topologies in {datetime.now() - start} seconds")
```

```{code-cell}
:tags: [hide-input]
ntips = len(sample_populations)
styles = ".sample text.lab {baseline-shift: super; font-size: 0.7em;}"
node_labels = {}

for p in range(ntips):
    name = ts.population(sample_populations[p]).metadata["id"]
    node_labels[p] = names[name]
    styles += f".n{p}>.sym {{fill: {colours[name]} }}"

total = sum(topology_span.values())
for rank, weight in topology_span.items():
    label = f"{weight/total *100:.1f}% of genome"
    embedded_tree = tskit.Tree.unrank(ntips, rank)
    display(embedded_tree.draw_svg(size=(160, 150), style="".join(styles), node_labels=node_labels, x_label=label))
```

Perhaps unsurprisingly, the most common topology is the one that groups the non-African
populations together (although there are many trees of the other two topologies,
mostly reflecting genetic divergence prior to the emergence of humans out of Africa).

For an example with real data, see {ref}`sec_popgen_topological`
in the {ref}`sec_intro_popgen` tutorial.