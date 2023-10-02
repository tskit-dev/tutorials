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

(sec_simulation_overview)=

# Tree sequences and simulation

**Yan Wong, Georgia Tsambos, and Peter Ralph**

Simulations are important in population genetics for many reasons:

::::{margin}
:::{todo}
Add links to papers that illustrate each of the following points
:::
::::

Exploration
: Simulations allow us to explore the influence of complex historical scenarios on
observed patterns of genetic variation and inheritance.

Benchmarking and evaluating methodologies
: To assess the accuracy of inferential methods, we need test datasets for which the
true values of important parameters are known.

Model training
: Some methods for ancestry inference are trained on simulated data (eg. Approximate
Bayesian Computation). This is especially important in studies of complex demographies,
where there are many potential parameters and models, making it impractical to specify
likelihood functions.

Compare to expectations
: It is often useful to compare data to what is expected under a simpler situation
(e.g. for use as a null model). For instance,  comparison to *neutral* simulations
can be used to identify regions subject to selection.

There are two major forms of population genetic simulation: **forwards-time**
and **backwards-time**. In general, forwards-time simulation is detailed and more
realistic, while backwards-time simulation is fast and efficient.

More specifically, apart from a
{ref}`few exceptions <msprime:sec_ancestry_models_selective_sweeps>`,
backwards-time simulations are primarily focused on neutral simulations, while
forward simulation is better suited to complex simulations, including those involving
selection and continuous space.

## Advantages of tree sequences

Some forwards-time ([SLiM](http://messerlab.org/slim/),
[fwdpy](http://molpopgen.github.io/fwdpy/)) and backwards-time
([msprime](https://tskit.dev/msprime)) simulators have a built-in capacity to output
tree sequences. This can have several benefits:

1. Neutral mutations, which often account for the majority of genetic variation, do not
    need to be tracked during the simulation, but can be added afterwards. See
    "{ref}`sec_tskit_no_mutations`".
2. Tree sequences can be used as an interchange format to combine backwards and
    forwards simulations, allowing you to take advantage of the advantages of both
    approaches. This is detailed in {ref}`sec_completing_forwards_simulations`.

## Some tips on simulation

Even with fast modern software, simulating full genome sequences of entire populations
can take some time. If you are finding your simulations too slow, it is worth
benchmarking them by running on a range of shorter chromosomes or sample sizes, then 
extrapolating to figure out how long the simulations you actually want to run would take.

:::{todo}
Add an example with a matplotlib fitted curve for some msprime simulations with
e.g. a high recombination rate.

Collecting data from simulations that take minutes to a few hours and looking at
the msprime paper for suggestions of what curve to fit to the data should give you
good predictions. See [issue #104](https://github.com/tskit-dev/tutorials/issues/104)
:::
