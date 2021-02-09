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

(sec_what_is)=

# What is a tree sequence?

A *succinct tree sequence*, or tree sequence for short, represents the relationships
between a set of DNA sequences. Tree sequences can be used to store genetic data
efficiently, and enable powerful analysis of millions of whole genomes at a time.
They can be created by [simulation](https://tskit.dev/software/#simulate) or by
[inferring relationships from genetic variation](https://tskit.dev/software/#infer).

As the name suggests, the simplest way to think about a tree sequences is as a set of
"local trees" - i.e. trees located at different points along the chromosome. Here's an
example based on ten genomes, $\mathrm{a}$ to $\mathrm{j}$, spanning a short 1kb
chromosome.

```{code-cell} ipython3
import string
import msprime
from IPython.display import SVG
```

```{code-cell} ipython3
seed = 3096  # chosen to create a nice example simulation
ts = msprime.sim_ancestry(5, population_size=1e4, sequence_length=1000,
    recombination_rate=1e-8, random_seed=seed)
# Extra code to label and order the tips alphabetically rather than numerically
genome_order = [n for n in ts.first().nodes(order="minlex_postorder") if ts.node(n).is_sample()]
labels = {n:string.ascii_lowercase[i] for i, n in enumerate(genome_order)}
style1 = ".node:not(.sample) > .sym {visibility: hidden;}"  # hide internal tree nodes
sz = (800, 250)  # size of the plot, slightly larger than the default

SVG(ts.draw_svg(size=sz, node_labels=labels, style=style1))
```

For almost three quarters of the chromosome, from the
start up to position 715, the relationships between the ten genomes are shown by
the first tree. The second tree shows the relationships between positions 715 and 932,
and the third from position 932 until the end.

These trees describe the full ancestry, or *genetic genealogy* of our 10 genomes. As
we shall see, this genealogy is [useful in many ways](sec_what_is_ancestry), and forms
the basis for a powerful [framework for statistical calculation](sec_what_is_analysis).
But first we shall see how the ancestry can be used as an efficient
*evolutionary encoding* for DNA sequences.

(sec_what_is_dna_data)=

## An efficient encoding of DNA data

A tree sequence can be used to describe patterns of genetic variation by combining the
trees with a knowledge of where *mutations* occur on their branches. Here's how that
might look in our simple example:

```{code-cell} ipython3
seed = 3
mutated_ts = msprime.sim_mutations(ts, rate=1e-7, random_seed=seed)

mut_labels = {}  # An array of labels for the mutations, listing position & allele change
l = "{:g} ({}→{})"
for mut in mutated_ts.mutations():  # This entire loop is just to make pretty labels
    site = mutated_ts.site(mut.site)
    older_mut = mut.parent >= 0  # is there an older mutation at the same position?
    prev = mutated_ts.mutation(mut.parent).derived_state if older_mut else site.ancestral_state
    mut_labels[mut.id] = l.format(site.position, prev, mut.derived_state)

SVG(mutated_ts.draw_svg(size=sz, node_labels=labels, mutation_labels=mut_labels, style=style1))
```

The trees tell us that, for example, the final mutation (at position 980) is inherited
by genomes $\mathrm{a}$ to $\mathrm{i}$. These genomes must have a *G* at that position,
compared to the original value of *C*, which only genome $\mathrm{j}$ has retained. In
other words, the combination of ancestry-plus-mutations fully defines the genetic
variation at all 11 variable sites. The variant at each of these 11 sites can be easily
extracted from the tree sequence:

```{code-cell} ipython3
haplotypes = mutated_ts.haplotypes()
print("\n".join(sorted([f"Genome {labels[i]}: {h}" for i, h in enumerate(haplotypes)])))
```

This is only a tiny example; real-world use cases may involve many millions of genomes
with hundreds of thousands of trees, spanning chromosomes of hundreds of megabases
in length. Tree sequences scale efficiently to such huge datasets because of one key
property: adjacent trees along the genome are highly correlated, that is, they
*share structure*. In our example, this becomes evident if we highlight the branches,
or "edges" in tree sequence terminology, that remain unchanged between the first
and the second tree (more visualization possibilities are detailed in the
[visualization tutorial](sec_tskit_viz)).

```{code-cell} ipython3
kept_edges = [e for e in ts.edges() if e.left==0 and e.right>ts.breakpoints(True)[1]]
style3 = (
    ",".join(f"#svg1 .tree:not(.t2) .node.a{e.parent}.n{e.child} > .edge" for e in kept_edges)
    + "{stroke:#00DD00; stroke-width: 2px}"
    + style1)
SVG(ts.draw_svg(size=sz, root_svg_attributes={'id':'svg1'}, node_labels=labels, style=style3))
```

Another way to think about shared structure is to notice that the second tree can be
formed by a simple rearrangement of the first tree. This can be done by simply switching
the centre group of five genomes, labelled $\mathrm{d}$ to $\mathrm{h}$, next to
$\mathrm{a}+\mathrm{b}+\mathrm{c}$. Similarly, the third tree just involves a single
adjustment: the movement of genome $\mathrm{i}$ away from being the closest relative of
$\mathrm{j}$. These sort of small rearrangements are typical of how genetic relationships
change along chromosomes, in both simulated and real datasets. Here's the take-home
message:

```{epigraph}
Tree sequences are efficient because they don't store each tree separately
```

More specifically, if adjacent trees share the same branch ("edge"), it need only be
stored once. As the trees get larger and larger, the efficency
gains become significant. For instance, if we simulate the evolution of a 10 billion
human chromosomes - roughly every living human - and store the resulting DNA sequences in
tree sequence format, it take 20,000 times less space than using the conventional
(so-called VCF) format and is about XXX times faster to process.

(sec_what_is_ancestry)=

## A record of genetic ancestry

Often, it turns out that what we are interested in is not the DNA sequence data, but the
genealogy itself. Knowledge of the ancestry can be used, for instance, to determine the
origin and age of variants under selection, to capture the spatial structure of
populations, or to uncover the effects of hybridization and admixture in the past.
Tree sequences are a particularly powerful representation because of how
closely they are linked to the underlying biological processes that generated the genomes
in the first place.

For example, each branch point in one of the trees above represents a most recent common
ancestor (MRCA), in other words a genome which existed at a specific time in the past.
It is helpful to distinguish these *ancestral* genomes from the *sampled* genomes
($\mathrm{a}$ to $\mathrm{j}$) which we have measured more directly. We can
indicate this in our trees by adding the MRCA genomes as circular nodes, rather than the
squares we have used for sampled genomes. 

```{note}
For clarity in these examples, we have been relabelling the sample genomes as
$\mathrm{a}$ to $\mathrm{j}$, and we will continue to do so. However, in standard tree
sequences, all the genomes (which are referred to as "nodes"), including the samples,
are numbered sequentially from 0. The ancestral nodes are in this example are therefore
labelled from 10 upwards.
```

```{code-cell} ipython3
all_labels = {n:n for n in range(ts.num_nodes)}
all_labels.update(labels)  # For consistency, keep letters for the sample nodes
tree_nodes = [set(tree.nodes()) for tree in ts.trees()]
kept_nodes = set.intersection(*tree_nodes)
changed_nodes = set(range(ts.num_nodes)) - kept_nodes
style2 = ",".join(f"#svg2 .n{n} > .sym" for n in changed_nodes) + "{fill:red}"
style2 += "#svg2 .node > .sym {visibility: visible;}"  # force-show all nodes: not normally needed
SVG(ts.draw_svg(size=sz, root_svg_attributes={'id':'svg2'}, node_labels=all_labels, style=style2))
```

Here we have also coloured those nodes (genomes) that only exist in some of the trees.
With longer chromosomes or higher recombination rates, essentially all the ancestral
nodes will fall into this category. Biologically, this means that most ancestral nodes
represent only a fraction of the whole ancestral genome that must have existed in the past.

The change from one tree to another is also biologically meaningful. It indicates that
one or more recombination events occured at this genomic location in the past. Note,
however, that for efficiency reasons and more, neither the recombination event itself
nor the branches on which it occurred are usually present in a tree sequence, although
it is possible to incorporate them via simulation (see the ARG tutorial).

(sec_what_is_analysis)=

## An efficient analysis framework

One of the most powerful aspects of the tree sequence concept is that all statistical
summaries of genetic variation data can be thought of in terms of trees (cite
Peter's paper).

Very brief intro to incremental algorithms

Very brief intro to combinatoric (topological) stuff

## How is a tree sequence stored

Under the hood, a tree sequence simply consists of a set of tables. Etc etc. Possibly a
picture of the edge spans?


## Why does it work?

Stuff here about philosophy and SPRs

Point out the similarity between a tree sequence and an ARG.
