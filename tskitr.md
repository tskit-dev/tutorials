---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.9.1
kernelspec:
  display_name: R
  language: R
  name: ir
---

# Tskit and R




```{code-cell}
library(reticulate)
msprime <- import("msprime")

ts <- msprime$sim_ancestry(5, random_seed=42)
ts
```

