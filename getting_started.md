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

(sec_tskit_getting_started)=
# Getting started with tskit

You've run some simulations or inference methods, and you now have a 
{class}`tskit.TreeSequence` object; what now? This tutorial is aimed 
users who are new to {program}`tskit` and would like to get some
basic tasks completed. We'll look at four fundamental things you might
need to do, and provide some pointers to where you can learn more.

First, let's simulate an example tree sequence using [msprime](https://tskit.dev/msprime):


```{code-cell} 
import msprime
 
ts = msprime.sim_ancestry(
    20,
    population_size=10_000,
    sequence_length=1_000_000,
    recombination_rate=1e-8)   
ts
```

## Processing trees

## Processing mutations

## Compute statistics

## Exporting data
