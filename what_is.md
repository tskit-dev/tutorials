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

```{code-cell} ipython3
:tags: [remove-cell]
import msprime

def whatis_example():
    seed = 3096  # chosen to create a nice example simulation
    ts = msprime.sim_ancestry(5, population_size=1e1, sequence_length=1000,
        recombination_rate=1e-5, random_seed=seed)
    # Mutate
    seed = 244  # a simpler example uses 214
    mutated_ts = msprime.sim_mutations(ts, rate=5e-5, random_seed=seed)
    
    mutated_ts.dump("data/whatis_example.trees")
    

def create_notebook_data():
    whatis_example()

# create_notebook_data()  # uncomment to recreate the tree seqs used in this notebook
```

(sec_what_is)=

# What is a tree sequence?

A *succinct tree sequence*, or "tree sequence" for short, represents the evolutionary
relationships between a set of DNA sequences. Tree sequences are based on fundamental
biological principles of inheritance, DNA duplication, and recombination; they can be
created by [simulation](https://tskit.dev/software/#simulate) or by
[inferring relationships from empirical DNA data](https://tskit.dev/software/#infer).

Tree sequences provide an efficient way of storing
[genetic variation](https://en.wikipedia.org/wiki/Genetic_variation) data, and can
power analyses of millions of whole [genomes](https://en.wikipedia.org/wiki/Genome).
Plots (a) and (b) summarize results presented
[further](plot_storing_everyone) [down](plot_incremental_calculation) this tutorial.

```{code-cell} ipython3
:"tags": ["hide-input"]
from IPython.display import SVG
import matplotlib_inline
import matplotlib.pyplot as plt
import numpy as np
%matplotlib inline
matplotlib_inline.backend_inline.set_matplotlib_formats('svg')

data1 = np.genfromtxt("data/storing_everyone.csv", delimiter=",", usecols=np.arange(1,12), names=True)
data2 = np.genfromtxt("data/benchmarks_without_copy_longer_genome.txt", encoding=None, names=True, dtype=None)
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(16, 4.5))
fig.subplots_adjust(wspace=0.5, left=0, right=1)
keep = data1['sample_size'] <= 1e6
x, y = data1['sample_size'][keep], data1['tsk_fit'][keep]/data1['vcf_fit'][keep]
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.loglog(x, y, c="C0", linewidth=4)
ax1.set_xlabel('# of 100Mb genomes', fontsize=18)
ax1.set_ylabel('Size of tree sequence\nfile (relative to VCF) ', fontsize=18)
ax1.tick_params(axis="both", labelsize=16)

txt = ax1.text(0.5, 1.3, "(a) Storing a million genomes as a tree sequence takes thousands of times less disk space",
    ha='center', va='top', transform=ax1.transAxes, wrap=True, size=24)
txt._get_wrap_line_width = lambda: 600

ts_time = {n: t for s, n, t in data2[['toolkit','nsam','seconds']] if s == 'tskit'}
libseq_time = {n: t for s, n, t in data2[['toolkit','nsam','seconds']] if s == 'libseq'}
x = np.unique(list(ts_time.keys()) + list(libseq_time.keys()))
y = np.array([libseq_time[time]/ts_time[time] for time in x])
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.loglog(x, y, linewidth=4)
ax2.set_xlabel("# of genomes", fontsize=18)
ax2.set_ylabel("Tajima's D calculations per\nunit time (relative to libseq)", fontsize=18)
ax2.tick_params(axis="both", labelsize=16)
txt = ax2.text(0.5, 1.3, "(b) Genetic calculations on millions of genomes can be sped up by many orders of magnitude",
    ha='center', va='top', transform=ax2.transAxes, wrap=True, size=24
)    
txt._get_wrap_line_width = lambda: 600
plt.show()
```

As the name suggests, the simplest way to think
about a tree sequence is as a sequence of "local trees" --- i.e. trees located at
different points along the [chromosome](https://en.wikipedia.org/wiki/Chromosome).
Here's a tiny example based on ten genomes, $\mathrm{a}$ to $\mathrm{j}$, spanning
a short 1000 letter chromosome.

```{code-cell} ipython3
:"tags": ["hide-input"]
import string
import tskit
from IPython.display import SVG

mutated_ts = tskit.load("data/whatis_example.trees")
ts = mutated_ts.delete_sites(list(range(mutated_ts.num_sites)))
# Extra code to label and order the tips alphabetically rather than numerically
labels = {i: string.ascii_lowercase[i] for i in range(ts.num_nodes)}
genome_order = [n for n in ts.first().nodes(order="minlex_postorder") if ts.node(n).is_sample()]
labels.update({n: labels[i] for i, n in enumerate(genome_order)})
style1 = (
    ".node:not(.sample) > .sym, .node:not(.sample) > .lab {visibility: hidden;}"
    ".mut {font-size: 12px}")
sz = (800, 250)  # size of the plot, slightly larger than the default

SVG(ts.draw_svg(
    size=sz, node_labels=labels, style=style1, time_scale="log_time", y_label="Time ago",
    y_axis=True, y_ticks=[0, 2, 5, 10, 20, 50, 100]))
```

The tickmarks on the X axis and background shading indicates the genomic positions covered
by the trees. For almost three quarters of the chromosome, from the
start until position 715, the relationships between the ten genomes are shown by
the first tree. The second tree shows the relationships between positions 715 and 932,
and the third from position 932 to the end. We can say that the first tree spans 715 base
pairs, the second 217, and the third 68.

Multiple trees are needed because of
[genetic recombination](https://en.wikipedia.org/wiki/Genetic_recombination), which causes
different regions of the chromosome to have different histories. Together, the sequence
of trees describe the full genetic ancestry, or *genetic genealogy*, of our 10 genomes.

(sec_what_is_dna_data)=

## An efficient encoding of DNA data

A tree sequence can be used to describe patterns of genetic variation by combining the
trees with a knowledge of where *mutations* occur on their branches. Here's how that
might look in our simple example:

```{code-cell} ipython3
:"tags": ["hide-input"]

mut_labels = {}  # An array of labels for the mutations, listing position & allele change
l = "{:g} ({}→{})"
for mut in mutated_ts.mutations():  # This entire loop is just to make pretty labels
    site = mutated_ts.site(mut.site)
    older_mut = mut.parent >= 0  # is there an older mutation at the same position?
    prev = mutated_ts.mutation(mut.parent).derived_state if older_mut else site.ancestral_state
    mut_labels[mut.id] = l.format(site.position, prev, mut.derived_state)

SVG(mutated_ts.draw_svg(
    size=sz, style=style1, time_scale="log_time",
    node_labels=labels, mutation_labels=mut_labels))
```

There are now ten single nucleotide mutations in the tree sequence. They are shown on the
branches of the trees, and the positions of the ten variable sites associated with the
mutations are shown along the X axis.

The trees inform us that, for example, the final mutation (at position 986) is inherited
by genomes $\mathrm{a}$ to $\mathrm{i}$. These genomes must have a *G* at that position,
compared to the original value of *C*. In other words, once we know the ancestry, placing
a relatively small number of mutations is enough to explain all the observed genetic
variation. Here's the resulting "variant matrix":

```{code-cell} ipython3
:"tags": ["hide-input"]
haplotypes = ["   ".join(h) for h in mutated_ts.haplotypes()]
print("Position: " + " ".join(str(int(s.position)) for s in mutated_ts.sites()))
print("\n".join(sorted(
    [f"Genome {labels[i]}:  {h}" for i, h in zip(mutated_ts.samples(), haplotypes)])))
```

This approach scales effectively to millions of genomes, and to chromosomes of
hundreds of megabases in length. The ability to deal with huge datasets comes down to
one key feature of genomic data: adjacent trees along a chromosome are highly correlated,
that is, they *share structure*. In our example this becomes evident
if we highlight the branches ("edges" in tree sequence terminology) that remain
unchanged between the first and the second tree.

```{code-cell} ipython3
:"tags": ["hide-input"]
# Highlight certain edges in certain trees. Other visualization possibilities in tutorials/viz.html
kept_edges = [e for e in ts.edges() if e.left==0 and e.right>ts.breakpoints(True)[1]]
style3 = (
    ",".join(f"#svg1 .tree:not(.t2) .node.a{e.parent}.n{e.child} > .edge" for e in kept_edges)
    + "{stroke:#00DD00; stroke-width: 2px}"
    + style1)
sz = (500, 250)
SVG(ts.draw_svg(
    size=sz, x_lim=(0, 900), root_svg_attributes={'id':'svg1'}, time_scale="log_time", y_axis=True,
    y_label="Time ago", y_ticks=[0, 2, 5, 10, 20, 50, 100], node_labels=labels, style=style3))
```

% Another way to think about shared structure is to notice that the second tree can be
% formed by a simple rearrangement of the first tree. This can be done by simply switching
% the centre group of five genomes, labelled $\mathrm{d}$ to $\mathrm{h}$, next to
% $\mathrm{a}+\mathrm{b}+\mathrm{c}$. Similarly, the third tree just involves a single
% adjustment: the movement of genome $\mathrm{i}$ away from being the closest relative of
% $\mathrm{j}$. These sort of small rearrangements are typical of how genetic relationships
% change along chromosomes, in both simulated and real datasets.
%%% possible link here to a tutorial which talks about SPRs

A branch can be shared by many adjacent trees, but is stored as a single edge in the tree
sequence. For large datasets this is a great saving, because typically each tree-change
affects only a few branches at a time, regardless of the tree size.
Here's the take-home message:

```{epigraph}
Tree sequences are efficient because they don't store each tree separately
```

Below is an extension of the plot at the top of this page, showing predicted
file sizes when storing not just millions, but billions of human-like genomes:
enough to encompass every human on the planet. This demonstrates that the tree sequence
encoding leads to savings of many orders of magnitude, even when compared against
compressed versions of the standard VCF storage format (original published data
[here](https://www.nature.com/articles/s41588-019-0483-y/figures/1)). It's also worth
noting that the efficiency extends to processing time too: tree sequences are often
several orders of magnitude faster to process than other storage formats.

(plot_storing_everyone)=

```{code-cell} ipython3
:"tags": ["hide-input"]
x = data1['sample_size']
fig, ax1 = plt.subplots(1, figsize=(10, 4))
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

plt.loglog(x,  data1['vcf_fit'], c="C1", label="VCF", linewidth=2)
plt.loglog(x,  data1['vcfz_fit'], c="C1", label="compressed VCF", linewidth=2, linestyle=":")

plt.loglog(x, data1['tsk_fit'], c="C0", label="tree sequence", linewidth=2)
plt.loglog(x, data1['tskz_fit'], c="C0", label="compressed tree sequence", linewidth=2, linestyle=":")

plt.xlabel('Number of 100Mb genomes (log scale)', fontsize=12)
plt.ylabel('Space required (GB, log scale)', fontsize=12)
plt.text(16e9, 0.001, 'Size of\nentire\nhuman\npopulation', ha="center", va="bottom", size=14)
plt.annotate('', xy=(16e9, 0.0001), xytext=(16e9, 0.001), 
            arrowprops=dict(facecolor='black', shrink=0))
plt.legend()
plt.show()
```

(sec_what_is_ancestry)=

## A record of genetic ancestry

Often, we're not interested so much in the DNA sequence data as the genetic ancestry
itself (as discussed in [this summary](https://www.nature.com/articles/s41588-019-0492-x)).
In other words, the main consideration is the actual trees in a tree sequence, rather
than the distributions of mutations placed upon them (indeed in genetic simulations, it
{ref}`may not be necessary<sec_tskit_no_mutations>` to incorporate neutral mutations at all).
The trees can be used, for example, to determine the origin and age of alleles under
selection, to capture the spatial structure of populations, or to uncover the effects
of hybridization and admixture in the past.

```{todo}
Insert illustration of the above, e.g. use of branch length calculations rather than
variants using colours for different branch lengths, or possibly a simple view of a
tree sequence over geographical space.
```


A major benefit of "tree sequence thinking" is the close relationship between the
tree sequence and the underlying biological processes that produced the genetic sequences
in the first place. For example, each branch point (or "internal node") in one of the
trees above can be imagined as a genome which existed at a specific time in the past, and
which is a "most recent common ancestor" (MRCA) of the descendant genomes at that
position on the chromosome. We'll mark these extra "ancestral genomes" on our picture,
distinguishing them from the *sampled* genomes ($\mathrm{a}$ to $\mathrm{j}$) by using
circular symbols: 

```{code-cell} ipython3
:"tags": ["hide-input"]
style2 = "#svg2 .node > .sym {visibility: visible;}"  # force-show all nodes: not normally needed
SVG(mutated_ts.draw_svg(
    size=sz, root_svg_attributes={'id':'svg2'}, y_label="Time ago",
    time_scale="log_time", y_axis=True, y_ticks=[0, 2, 5, 10, 20, 50, 100],
    node_labels=labels, mutation_labels={}, style=style2))
```

Knowing the tree sequence means that we can easily deduce the ancestral genomes
$\mathrm{k}$ to $\mathrm{u}$, by looking at which mutations they have inherited.

```{code-cell} ipython3
:"tags": ["hide-input"]
tables = mutated_ts.dump_tables()
# Flip sample and nonsample flags, making the haplotypes() method print out nonsample nodes
s_flags = tables.nodes.flags[ts.samples()[0]]
no_flags = s_flags-s_flags
tables.nodes.flags = np.where(tables.nodes.flags & tskit.NODE_IS_SAMPLE, no_flags, s_flags)
ts_flipped = tables.tree_sequence()
haplotypes = ["   ".join(h) for h in ts_flipped.haplotypes(missing_data_character=" ")]
print(" " * ts_flipped.num_sites, " " * (ts_flipped.num_sites-4), "")
print(
    "||ANCESTRAL GENOMES||    Position:",
    " ".join(str(int(s.position)) for s in ts_flipped.sites()))
print(
    "\n".join(reversed(sorted([
        f"Genome {labels[i]} (time {ts.node(i).time:5.2f} in the past):  {h}"
        for i, h in zip(ts_flipped.samples(), haplotypes)]))))
```
You can see that some ancestors (particularly the older ones) are missing genomic regions,
because those parts of their genome have not been inherited by any of the sampled
genomes (i.e. that ancestral node is not part of the tree at that position in the sequence)

```{note}
For clarity in these examples, we have been using letters to label nodes, but normally
the nodes are referred to by number.
```

```{todo}
Mention ARGs in passing and link out to the ARG tutorial.
% Somewhere we should explain *why* trees change along the genome, and it
% would be good to mention ARGs in passing somewhere. We previously had too much
% detail, though:
%
% The change from one tree to another is also biologically meaningful. It indicates that
% one or more recombination events occured at this genomic location in the past. Note,
% however, that for efficiency reasons and more, neither the recombination event itself
% nor the branches on which it occurred are usually present in a tree sequence, although
% it is possible to incorporate them via simulation (see the ARG tutorial).
```

(sec_what_is_analysis)=

## A framework for efficient computation


```{todo}
Introduction: algorithms on trees are known to be efficient (phylogenetics). We
extend these to multiple correlated trees. Mention "dynamic programming" in passing.
```

Statistical measures of genetic variation can be thought of as a calculation combining
the local trees with the mutations on each branch (or, often preferably, the length of the
branches: see [this summary](https://www.genetics.org/content/genetics/215/3/779)).
Because a tree sequence is built on a set of small branch changes along the chromosome,
statistical calculations can often be updated incrementally as we
move along the genome, without having to perform the calculation *de novo* on each tree.
Using {program}`tskit`, the tree sequence toolkit, can result in speed-ups of many
orders of magnitude when perfoming calculations on large datasets, as in this example of
calculating [Tajima's D](https://en.wikipedia.org/wiki/Tajima%27s_D)
(from [here](https://www.genetics.org/content/215/3/779#F9)):

(plot_incremental_calculation)=
```{code-cell} ipython3
:"tags": ["hide-input"]
ts_time = np.array([[n,t] for s, n, t in data2[['toolkit','nsam','seconds']] if s == 'tskit'])
ska_time = np.array([[n, t] for s, n, t in data2[['toolkit','nsam','seconds']] if s == 'allel'])
libseq_time = np.array([[n, t] for s, n, t in data2[['toolkit','nsam','seconds']] if s == 'libseq'])
fig, ax1 = plt.subplots(1, figsize=(10, 5))
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.loglog(ska_time[:,0], ska_time[:,1], c="C3", linewidth=2, label="scikit-allel library")
ax1.loglog(libseq_time[:,0], libseq_time[:,1], c="C1", linewidth=2, label="libseq library")
ax1.loglog(ts_time[:,0], ts_time[:,1], c="C0", linewidth=2, label="tskit library")
ax1.set_ylabel("Time to calculate Tajima's D (secs/site)", fontsize=12)
ax1.set_xlabel("Number of sampled genomes", fontsize=12)
plt.legend()
plt.show()
```

```{todo}
Very brief discussion of efficient counting of topologies, i.e. the combinatorics module
```

Summary of this subsection:

```{epigraph}
Genetic calculations involve iterating over trees, which is highly efficient in tskit 
```


## Further reading

* How is a tree sequence stored: details in the
  [data structures](sec_data_structures) tutorial
* The offical {program}`tskit` [documentation](https://tskit.dev/tskit/docs)
  [data structures](sec_data_structures) tutorial
* The tree sequence philosophy. biological underpinnings and SPRs (to do)
