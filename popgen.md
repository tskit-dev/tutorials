---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

```{currentmodule} tskit
```

```{code-cell}
:tags: [remove-cell]
import urllib.request

import tqdm
import tskit
import tszip

class DownloadProgressBar(tqdm.tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def download(url, progress=True):
    with DownloadProgressBar(
        unit='B',
        unit_scale=True,
        miniters=1,
        desc=url.split('/')[-1],
        disable=not progress,
    ) as t:
        tmp_fn, _ = urllib.request.urlretrieve(url, reporthook=t.update_to)
        try:
            ts = tskit.load(tmp_fn)
        except tskit.FileFormatError:
            # could be a tsz file
            ts = tszip.decompress(tmp_fn)
        urllib.request.urlcleanup() # Remove tmp_fn
    return ts

def download_unified_genealogy():
    keep_span = [108_000_000, 110_000_000]  # cut down to this genome region
    keep_regions = {"EastAsia", "EAST_ASIA", "AFRICA", "Africa"}

    # Downloads 138 Mb of data - this may take a while
    tables = download(
        "https://zenodo.org/record/5512994/files/"
        "hgdp_tgp_sgdp_high_cov_ancients_chr2_q.dated.trees.tsz"
    ).dump_tables()
    tables.keep_intervals([keep_span])
    tables.populations.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.sites.metadata_schema = tskit.MetadataSchema.permissive_json()
    ts = tables.tree_sequence()
    ts = ts.simplify([
        u
        for u in ts.samples()
        if (
            ts.population(ts.node(u).population).metadata.get("region") in keep_regions
            or ts.population(ts.node(u).population).metadata.get("name") == "Denisovan"
        )
    ])
    tszip.compress(ts, "data/unified_genealogy_2q_108Mb-110Mb.tsz")

def create_notebook_data():
    download_unified_genealogy()

# create_notebook_data()  # uncomment to recreate the tree seqs used in this notebook
```

(sec_intro_popgen)=

# `Tskit` for population genetics

{ref}`Tskit<tskit:sec_introduction>`, the tree sequence toolkit, brings the power of
evolutionary trees to the field of population genetics. The
{ref}`succinct tree sequence<sec_what_is>` format
is designed to store DNA sequences jointly with their ancestral history (the
"genetic genealogy" or {ref}`ARG<sec_args>`). Storing population genetic data in this
form enables highly efficient computation and analysis.

The core `tskit` library provides methods for storing genetic data, a flexible
analysis framework, and APIs to build your own efficient population genetic algorithms.
Because of its speed and scalability, `tskit` is well-suited to interactive analysis of
large genomic datasets. 

## Population genetic simulation

Several simulation tools output tree sequences. Below we use the
standard library for population genetic simulation models
([stdpopsim](https://popsim-consortium.github.io/stdpopsim-docs/)) to generate a model of
*Homo sapiens*, in which African, Eurasian,
and Asian populations combine to generate a mixed American population. We can use the
[demesdraw](https://pypi.org/project/demesdraw/) package to plot a schematic of the
migrations and population size changes that define this model. 


```{code-cell}
import stdpopsim
import demesdraw
from matplotlib import pyplot as plt

species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model("AmericanAdmixture_4B18")

# Plot a schematic of the model
demesdraw.tubes(model.model.to_demes(), ax=plt.gca(), seed=1, log_time=True)
plt.show()
```

Genomic data in tree sequence format can be generated via the widely-used
[msprime](https://tskit.dev/software/msprime.html) simulator. Here we simulate 1
megabase of genome sequence at the start of human chromosome 1 under this model,
together with its evolutionary history. We generate 16 diploid genomes: 4 from each of
the populations in the model. The DNA sequences and their ancestry are stored in a
succinct tree sequence named `ts`:

```{code-cell}
contig = species.get_contig("chr1", mutation_rate=model.mutation_rate, right=1_000_000)
samples = {"AFR": 4, "EUR": 4, "ASIA": 4, "ADMIX": 4} # 16 diploid samples
engine = stdpopsim.get_engine("msprime")
ts = engine.simulate(model, contig, samples, seed=9).trim()  # trim to first 20kb simulated
print(f"Simulated a tree sequence of {ts.num_samples} haploid genomes:")
print(f"{ts.num_sites} variable sites over {ts.sequence_length} base pairs")
```

We can now inspect alleles and their frequencies at the variable sites we have simulated
along the genome:

```{code-cell}
for v in ts.variants():
    print(
        f"Variable site {v.site.id} at position {v.site.position} has allele frequencies",
        {state: f"{freq:.1%}" for state, freq in v.frequencies().items()}
    )
    if v.site.id > 4:
        print("...")
        break
```

Or we can efficiently grab the genotypes for each sampled genome

```{code-cell}
print("Sample ---> ", " ".join([f"{u:>2}" for p in ts.populations() for u in ts.samples(population=p.id)]))
print("Population |", "".join([f"{p.metadata['name']:^{3*(len(ts.samples(population=p.id)))-1}}|" for p in ts.populations()]))
print("__________ |", "".join(["_" * (3 * len(ts.samples(population=p.id)) - 1) + "|" for p in ts.populations()]))
print("  Position")
for v in ts.variants():
    print(f"{int(v.site.position):>10} | ", "  ".join(v.states()))
    if v.site.id >= 30: #  Only show the first 30 sites, for brevity
        break
```

It is also possible to grab the [haplotypes](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.haplotypes)
for specific samples (although this is slightly less efficient)

```{code-cell}
pop_id = 0
samples = ts.samples(population=pop_id)

for sample_id, haplotype in zip(samples, ts.haplotypes(samples=samples)):
    h = ".".join(list(haplotype))  # Add a dot between letters, to clarify they are not  adjacent
    print(f"Sample {sample_id:<2} ({ts.population(pop_id).metadata['name']:^5}): {h}")
```

You can easily obtain the
{meth}`TreeSequence.allele_frequency_spectrum` for the entire region (or for
{ref}`windowed regions<sec_tskit_getting_started_compute_statistics_windowing>`) 
directly from the tree sequence (i.e. without needing to reconstruct genotypes)

```{code-cell}
afs = ts.allele_frequency_spectrum()
plt.bar(range(ts.num_samples + 1), afs)
plt.title("Allele frequency spectrum")
plt.show()
```

Similarly `tskit` allows fast and easy
{ref}`calculation of statistics<sec_tutorial_stats>` along the genome. Here is
a plot of windowed $F_{st}$ between Africans and admixed Americans over this
region of chromosome:

```{code-cell}
# Define the samples between which Fst will be calculated
pop_id = {p.metadata["name"]: p.id for p in ts.populations()}
sample_sets=[ts.samples(pop_id["AFR"]), ts.samples(pop_id["ADMIX"])]

# Do the windowed calculation, using windows of 10 kilobases
windows = list(range(0, int(ts.sequence_length + 1), 10_000))
F_st = ts.Fst(sample_sets, windows=windows)

# Plot
plt.stairs(F_st, windows, baseline=None)
plt.ylabel("AFR-ADMIX Fst")
plt.xlabel("Genome position")
plt.show()
```

Extracting the genetic tree at a specific genomic location is easy using `tskit`, which
also provides methods to {ref}`plot<sec_tskit_viz>` these trees. Here we
grab the tree at position 10kb, and colour the different populations by
grab the tree at position 10kb, and colour the samples according to their population,
as described in the {ref}`viz tutorial<sec_tskit_viz_styling>`:

```{code-cell}
tree = ts.at(10_000)

colours = dict(AFR="yellow", EUR="cyan", ASIA="green", ADMIX="red")
styles = [
    f".leaf.p{pop.id} > .sym {{fill: {colours[pop.metadata['name']]}}}"
    for pop in ts.populations()
]

styles += [ # rotate the population labels, etc
    ".leaf > .lab {text-anchor: start; transform: rotate(90deg) translate(6px)}",
    ".leaf > .sym {stroke: black}"
]

labels = { # Label samples by population
    u: ts.population(ts.node(u).population).metadata["name"].capitalize()
    for u in ts.samples()
}

tree.draw_svg(
    size=(800, 500),
    canvas_size=(800, 520),
    node_labels=labels,
    style="".join(styles),
    y_axis=True,
    y_ticks=range(0, 30_000, 10_000)
)
```

Or we can plot a principal components analysis of the genome, which should reflect
geographical distinctiveness:

```{code-cell}
from matplotlib.patches import Patch

# Run the Principal Components Analysis (PCA)
pca_obj = ts.pca(num_components=2)

# Plot the PCA "factors"
col_list = [colours[pop.metadata["name"]] for pop in ts.populations()]
sample_pop_ids = ts.nodes_population[ts.samples()]
plt.scatter(*pca_obj.factors.T, c=[col_list[p] for p in sample_pop_ids], edgecolors= "black")
plt.xlabel("PCA 1")
plt.ylabel("PCA 2")
plt.legend(handles=[
    Patch(color=col_list[pop.id], label=pop.metadata["name"]) for pop in ts.populations()
]);
```

## Population genetic inference

If, instead of simulations, you want to analyse existing genomic data (for example
stored in a VCF file), you will need to infer a tree sequence from it, using e.g.
[tsinfer](https://tskit.dev/tsinfer/docs/stable/). Here we load an illustrative portion
of an [inferred tree sequence](https://zenodo.org/record/5512994)
based on about 7500 public human genomes, including genomes from the
[Thousand Genomes Project](https://www.internationalgenome.org/data-portal/data-collection/grch38) and
[Human Genome Diversity Project](https://www.internationalgenome.org/data-portal/data-collection/hgdp).
The genomic region encoded in this tree sequence has been cut down to
span positions 108Mb-110Mb of human chromosome 2, which spans the
[EDAR](https://en.wikipedia.org/wiki/Ectodysplasin_A_receptor) gene.

Note that we are using {func}`tszip:tszip.load` to load the file, as this
utility can also read and write compressed tree sequences in `.tsz` format.

```{code-cell}
import tszip
ts = tszip.load("data/unified_genealogy_2q_108Mb-110Mb.tsz")

# The ts encompasses a region on chr 2 with an interesting SNP (rs3827760) in the EDAR gene
edar_gene_bounds = [108_894_471, 108_989_220]  # In Mb from the start of chromosome 2
focal_variant = [v for v in ts.variants() if v.site.metadata.get("ID") == "rs3827760"].pop()
print("An interesting SNP within the EDAR gene:")
focal_variant
```

For simplicity, this tree sequence has been {ref}`simplified<sec_simplification>` to
include only those samples from the African and East Asian regions. These belong to a
number of populations. The population information, as well as information describing the
variable sites, is stored in tree sequence {ref}`metadata<sec_tutorial_metadata>`:

```{code-cell}
import pandas as pd

print(ts.num_populations, "populations defined in the tree sequence:")

pop_names_regions = [
    [p.metadata.get("name"), p.metadata.get("region"), len(ts.samples(population=p.id))]
    for p in ts.populations()
]
with pd.option_context('display.max_rows', 100):
    display(pd.DataFrame(pop_names_regions, columns=["name", "region", "# genomes"]))
```

You can see that there are multiple African and East asian populations, grouped by
region. Here we collect two lists of IDs for the sample
{ref}`nodes<sec_terminology_nodes>` from the African region and from the East asian
region:

```{code-cell}

sample_lists = {}
for n, rgns in {"Africa": {'AFRICA', 'Africa'}, "East asia": {'EAST_ASIA', 'EastAsia'}}.items():
    pop_ids = [p.id for p in ts.populations() if p.metadata.get("region") in rgns]
    sample_lists[n] = [u for p in pop_ids for u in ts.samples(population=p)]
```


With these lists we can calculate different windowed statistics
(here {meth}`genetic diversity<TreeSequence.diversity>` and
{meth}`Tajima's D<TreeSequence.Tajimas_D>`) within each of these regions:

```{code-cell}
edar_ts = ts.trim()  # remove regions with no data (changes the coordinate system)
windows = list(range(0, int(edar_ts.sequence_length)+1, 10_000))
data = {
    "Genetic diversity": {
        region: edar_ts.diversity(samples, windows=windows)
        for region, samples in sample_lists.items()
    },
    "Tajima's D": {
        region: edar_ts.Tajimas_D(samples, windows=windows)
        for region, samples in sample_lists.items()
    },  
}

# Plot the `data`
fig, axes = plt.subplots(ncols=2, figsize=(15, 3))
start = ts.edges_left.min()  # the empty amount at the start of the tree sequence

for (title, plot_data), ax in zip(data.items(), axes):
    ax.set_title(title)
    ax.axvspan(edar_gene_bounds[0], edar_gene_bounds[1], color="lightgray")
    ax.axvline(focal_variant.site.position, ls=":")
    for label, stat in plot_data.items():
        ax.stairs(stat, windows+start, baseline=None, label=label)
    ax.text(edar_gene_bounds[0], 0, "EDAR")
    ax.legend()
plt.show()
```

Other population genetic libraries such as
[scikit-allel](https://scikit-allel.readthedocs.io/en/stable/) (which is
{ref}`interoperable<sec_tskit_getting_started_exporting_data_allel>` with `tskit`)
could also have been used to produce the plot above. In this case, the advantage of
using tree sequences is simply that they allow these sorts of analysis to
{ref}`scale<plot_incremental_calculation>` to datasets of millions of whole genomes.

(sec_popgen_topological)=

### Topological analysis

As this inferred tree sequence stores (an estimate of) the underlying
genealogy, we can also derive statistics based on genealogical relationships. You
may have noticed that this tree sequence also contains a sample genome based on an ancient
genome, a [Denisovan](https://en.wikipedia.org/wiki/Denisovan) individual. We'll first
simplify the tree sequence to focus on only the Denisovan plus
a common East Asian and a common African population: 

```{code-cell}
# Focus on Han, San, and Denisovan
focal = {
    "Han": ts.samples(population=6),
    "San": ts.samples(population=17),
    "Denisovan": ts.samples(population=66),
}

for name, nodes in focal.items():  # Sanity check that we got the right IDs
    assert ts.population(ts.node(nodes[0]).population).metadata["name"] == name

# Simplify to just those samples ...
all_focal_samples = [u for samples in focal.values() for u in samples]
simplified_ts = ts.simplify(all_focal_samples, filter_sites=False)

# ... and find the tree around the rs3827760 SNP
focal_site = simplified_ts.site(focal_variant.site.id)
tree = simplified_ts.at(focal_site.position)
```

With this smaller number of samples, we can easily plot the tree
at the "rs3827760" SNP:

```{code-cell}
:"tags": ["hide-input"]
# Make some nice labels, colours, and legend
mutation_labels = {m.id: focal_site.metadata.get("ID") for m in focal_site.mutations}
colours = dict(San="yellow", Han="green", Denisovan="magenta")
styles = [
    f".leaf.p{pop.id} > .sym {{fill: {colours[pop.metadata['name']]}; stroke: grey}}"
    for pop in simplified_ts.populations()
]
legend = '<rect width="125" height="75" x="100" y="30" fill="transparent" stroke="grey" />'
legend += '<text x="120" y="45" font-weight="bold">Populations</text>'
# Create the legend lines, one for each population. Setting classes that match those
# used for normal nodes means that styled colours are auto automatically picked-up.
legend += "".join([
    f'<g transform="translate(105, {60 + 15*p.id})" class="leaf p{p.id}">'  # an SVG group
    f'<rect width="6" height="6" class="sym" />'  # Square symbol
    f'<text x="10" y="7">{p.metadata["name"]}' # Label
    f'{(" (" + p.metadata["region"].replace("_", " ").title() + ")") if "region" in p.metadata else ""}</text></g>'  
    for p in simplified_ts.populations()
])

tree.draw_svg(
    size=(1000, 400),
    style="".join(styles),
    node_labels={},
    mutation_labels=mutation_labels,
    preamble=legend,
    title=f"Tree of human chromosome 2 at position {int(focal_variant.site.position)}",
    y_axis=True,
    y_ticks=range(0, 50_000, 10_000),
)
```

You can see that the pair of magenta Denisovan genomes in this region tend to be
more closely associated with the East Asian genomes. We can assess that by counting
all the 3-tip topologies in the tree that contain one genome from each population:

```{code-cell}
topology_counter = tree.count_topologies()
embedded_topologies = topology_counter[range(simplified_ts.num_populations)]
```

```{code-cell}
:"tags": ["hide-input"]
# All the following code is simply to plot the embedded_topologies nicely
all_trees = list(tskit.all_trees(simplified_ts.num_populations))
last = len(all_trees) - 1
svgs = ""
style = "".join(styles) + ".sample text.lab {baseline-shift: super; font-size: 0.7em;}"
style = style.replace(".leaf.p", ".leaf.n")  # Hack to map node IDs to population colours
params = {
    "size": (160, 150),
    "node_labels": {pop.id: pop.metadata["name"] for pop in simplified_ts.populations()}
}
for i, t in enumerate(all_trees):
    rank = t.rank()
    count = embedded_topologies[rank]
    params["title"] = f"{count} trees"
    if i != last:
        svgs += t.draw_svg(root_svg_attributes={'x': (last - i) * 150}, **params)
    else:
        # Plot the last svg and stack the previous ones to the right
        display(t.draw_svg(preamble=svgs, canvas_size=(1000, 150), style=style, **params))
```


See {ref}`sec_counting_topologies` for an introduction to topological methods in
`tskit`.

## Further information

This brief introduction is meant as a simple taster. Many other efficient population
genetic {ref}`analyses<sec_analysing_tree_sequences>` are possible when you have
genomic data stored as a tree sequence.

The rest of the {ref}`tutorials<sec_intro>` contain a large number of examples which
are relevant to population genetic analysis and research. You can also visit the
[learning section](https://tskit.dev/learn/) of the [tskit website](https://tskit.dev/).
