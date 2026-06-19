---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# PCA, on branches and SNPs

+++

Principal Component Analysis (PCA) is commonly used for exploring population structure in genetic datasets, where it is usually computed from SNP genotyped data. In the context of ARGs, it is also possible to perform branch PCA as implemented in `tskit`. This does not use variant data. (Of course, it may indirectly rely on variant data if the ARG was inferred from data.)

In this tutorial, we demonstrate both approaches. We will apply these to haplotypes and diploid genotype data.

The documentation of `tskit.TreeSequence.pca` can be found [here](https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.pca).

+++

:::{note}
Usually, PCA is carried out on a diploid genotype matrix (individuals in rows, loci in columns) with values 0, 1, and 2. PCA can then be achieved through singular value decomposition (SVD) of the column-centred genotype matrix. This results in a matrix of principal component (PC) scores, which are linear combinations of the genotype columns. The PC scores are ordered, decreasingly, by the amount of variation from the original data they account for.
:::

+++

First, we'll simulate an ARG with population structure:

```{code-cell} ipython3
# load required libraries
import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
```

```{code-cell} ipython3
# set a mutation rate
mu = 1e-8
# number of sub-populations/'islands'
nPop = 5
# number of diploids sampled from each sub-population
nSamp = 10
# number of haplotypes sampled
nHap = 2* nSamp
# per-island effective population size
nn = 1e4 
# migration rate (per individual and island pair)
migRate = 1e-5
```

Simulate an ARG using an island model demography. There are five islands, each with a population size of 10,000. Pairwise migration rates are $10^{-5}$.

```{code-cell} ipython3
# Island model demography, 5 islands connected by low gene flow
dmg = msprime.Demography.island_model([nn] * nPop, migration_rate=migRate)
```

```{code-cell} ipython3
# Simulate ARG
ts = msprime.sim_ancestry(samples={i: nSamp for i in range(nPop)},
                      demography=dmg,
                      random_seed=1234,
                      sequence_length=1e6,
                      recombination_rate=1e-8)
ts
```

The same ARG, but with mutations added.

```{code-cell} ipython3
# Add mutations
tsm = msprime.sim_mutations(ts, rate=mu, random_seed=1234)
tsm
```

The migration rates between the islands are quite low. This should lead to considerable genetic differentiation. Let us compute pairwise $F_{ST}$:

```{code-cell} ipython3
# Considerable pairwise Fst between the 'islands'
fstmat = np.zeros([nPop,nPop])
for i in range(nPop-1):
    for j in range(i+1,nPop):
        fstmat[i,j] = tsm.Fst([range(i*nHap,(i+1)*nHap), range((i+1)*nHap,(i+2)*nHap)])
fstmat
```

## Branch PCA (tskit)
To demstrate that branch PCA works without variant data, we run it on the ARG without mutations, `ts`.

```{code-cell} ipython3
# haplotypes, each sample haplotype is ues by default
hapBranchPca=ts.pca(num_components=10)
```

```{code-cell} ipython3
# genotypes, all individuals are specified
dipBranchPca=ts.pca(num_components=10, individuals=range(5*nSamp))
```

## PCA 'by hand'
To compute a traditional SNP PCA, we start by extracting the haploid 'genotypes' from the ARG. We then make use of the `TreeSequence` object's `individuals_nodes` property (an array) to select each individual's two haplotypes and to add them to create individual diploid genotypes.

```{code-cell} ipython3
# obtain a haplotype matrix from the tree sequence with mutation; print its shape
# 100 haplotypes (= 10 individual samples * 5 islands * 2 haplotypes per individual)
# 13683 variant sites
htMat=tsm.genotype_matrix().transpose()
htMat.shape
```

```{code-cell} ipython3
# Add each individual's two haplotypes to generate individual genotypes
sample_ids_to_mat_index = np.full_like(tsm.samples(), tskit.NULL, shape=tsm.num_nodes)
sample_ids_to_mat_index[tsm.samples()] = np.arange(len(tsm.samples()))
gtMat = htMat[sample_ids_to_mat_index[ts.individuals_nodes]].sum(axis=1)
```

```{code-cell} ipython3
# Haplotype SVD (column-centred)
hapSvd = np.linalg.svd(htMat - htMat.mean(axis=0), full_matrices=False)
```

```{code-cell} ipython3
# Genotype SVD (column-centred)
dipSvd = np.linalg.svd(gtMat - gtMat.mean(axis=0), full_matrices=False)
```

## Plot for comparison
Note that PCA does not preserve the axis orientation. The plots in the panels below will show similar patterns but one or both axes may be flipped.

```{code-cell} ipython3
fig, axs = plt.subplots(2, 2)
plt.tight_layout()
axs[0, 0].scatter(hapSvd.U[:,0],
                  hapSvd.U[:,1],
                  c=np.repeat([1,2,3,4,5], [nHap] * nPop))
axs[0, 0].set_title('Haplotypes (sites)')

axs[0,0].set_ylabel("PC2")
axs[0,1].scatter(dipSvd.U[:,0],
                 dipSvd.U[:,1],
                 c=np.repeat([1,2,3,4,5], [nSamp] * nPop))
axs[0,1].set_title("Individuals (sites)")

# flipping the axes to make similarity clearer:
axs[1,0].scatter(hapBranchPca.factors[:,0],
                 hapBranchPca.factors[:,1],
                 c=np.repeat([1,2,3,4,5], [nHap] * nPop))
axs[1,0].set_title("Haplotypes (branches)")
axs[1,0].set_ylabel("PC2")
axs[1,0].set_xlabel("PC1")

axs[1,1].scatter(dipBranchPca.factors[:,0],
                 dipBranchPca.factors[:,1],
                 c=np.repeat([1,2,3,4,5], [nSamp] * nPop))
axs[1,1].set_title("Individuals (branches)")
axs[1,1].set_xlabel("PC1")

plt.show()
```

The plots on the left show one dot per haplotype. These have twice as many dots as the plots on the right, which show individuals. The colours indicate from which of the five islands a haplotype or individual was sampled. As expected with low geneflow, there is some grouping by island. Feel free to re-run with higher or lower values of `migRate` to see how the separations between the island samples changes.

+++

## Comparing variance components between branch and SNP PCA
Both `numpy.linalg.svd` and `tskit.TreeSequence.pca` return information about the amount of variation accounted for by each PC. These information are stored in the slots `S` (standard variation for SVD) and `eigenvalues` (variance for branch PCA). To make the two match, we need to multiply the eigenvalues by the mutation rate before taking the square root.

```{code-cell} ipython3
# square root of (branch eigenvalues multiplied by the mutation rate)
xx=np.sqrt(hapBranchPca.eigenvalues * mu)
# SVD S values
yy=hapSvd.S[:10]
```

We now fit a least-squares regression model to demonstrate the match between SVD standard variation and transformed eigenvalues.

```{code-cell} ipython3
res = linregress(xx, yy)
print(f"Intercept: {res.intercept:.4f}\n    Slope: {res.slope:.4f}\n      r^2: {res.rvalue**2:.4f}")
```

$r^2$ is close to 1. Let us visualise this. Each dot below shows a standard deviation value associate with one PC. The fact that they are well correlated suggests that both SNP and branch PCA yielded very similar results.

```{code-cell} ipython3
plt.scatter(xx, yy)
plt.xlabel("sqrt(Branch eigenvals)")
plt.ylabel("GT svd.S")
plt.plot(xx, res.intercept + res.slope*xx, 'r', label='fitted line')
plt.xlabel(r"Branches: $\sqrt{eigenvals * \mu}$") # use raw string to avoid error message about \s
plt.ylabel("SNPs: $S$")
plt.title("Variance components of SNP and branch PCA")
plt.grid()
plt.show()
```

## Time windows
Above we showed how variant and branch-based PCA are equivalent. But the ARG is a much richer data type than the genotype matrix. ARGs contain information about the historic relationships between the samples (possibly blurred by a inference step). Branch PCA allows one to specify a time window over which the PCA is to be computed, something that cannot be done for SNP PCA. Next, we compute PCA in time slices with breaks 0, 10, 100, 1000, 10,000, 100,000, 100,0000, 1,000,000, and 10,000,000. The results are stored in a list.

```{code-cell} ipython3
pctime=[tsm.pca(num_components=10, time_windows=[10**i, 10**(i+1)]) for i in range(8)]
```

Being of class `PCAResult`, the elements of the list have a `factors` property. This has a shape of (100,10). I.e., 10 PCs for 100 haplotypes.

```{code-cell} ipython3
pctime[0].factors.shape
```

```{code-cell} ipython3
for i in range(8):
    evecs = pctime[i].factors[:,:10]

    plt.scatter(evecs[:,0],
                evecs[:,1],
                c=np.repeat(range(5), 20))
    plt.title(f"Branch PCA, (window: {10**i} - {10**(i+1)} gens)")
    plt.show()
```

When selecting a very old window, each individual contributes to its own PC, causing most to be plotted at the origin (0,0). We can see this when inspecting the oldest window's PC scores, which are an identityt matrix. All haplotypes below the first two have 0 entries for the first two PC scores (the two left-most columns).

```{code-cell} ipython3
pctime[7].factors[:20,:10]
```

## Empirical data
Here, we demonstrated using simulated data how SNPs and ARG branches lead to equivalent PCA results. For empirical data, the ancestral states of variant sites are not known a priori, which will in practice often lead to polarisation differences. That may affect the outcome of PCA.

**TODO:** Extend Tutorial to empirical data.
