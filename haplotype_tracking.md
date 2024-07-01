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


(sec_tutorial_haplotypes)=

# Tracking inheritance of haplotypes

A tree sequence provides an encoding of how segments of genome are inherited.
For some purposes, it is most helpful to iterate along the genome, looking
sequentially at each of the genealogical trees implied by this pattern of inheritance.
However, the data structure itself was not really designed for this purpose:
it naturally arose from the perspective of looking back through time, to see
how genomes were inherited from each other (in other words, the *coalescent* perspective).
This tutorial demonstrates how to use the information in the tree sequence to see
how these inherited segments of ancestry change as one moves through time.

To do this, it will be helpful define a simple class to represent a collection
of non-overlapping intervals. Each ancestral lineage will have an associated
collection of intervals, and as we move back in time we will update these.

```{code-cell} ipython3
class SegmentList:
    
    def __init__(self, segments=None):
        if segments is None:
            segments = []
        self.segments = segments
    
    def __str__(self):
        return self.segments.__str__()
    
    def __repr__(self):
        return self.segments.__repr__()
    
    def remove_segment(self, left, right):
        """
        Removes the intersection with [left, right).
        """
        removed = []
        new = []
        for a, b in self.segments:
            u, v = a, min(b, left)
            if u < v:
                new.append([u, v])
            u, v = max(a, left), min(b, right)
            if u < v:
                removed.append([u, v])
            u, v = max(a, right), b
            if u < v:
                new.append([u, v])
        self.segments = new
        return SegmentList(removed)
    
    def add_segment(self, left, right):
        """
        Updates to the union with [left, right).
        """
        new = []
        for a, b in self.segments:
            if max(a, left) < min(b, right):
                # overlaps
                left = min(left, a)
                right = max(right, b)
            else:
                new.append([a, b])
        new.append([left, right])
        self.segments = new
```

Here is a small tree sequence:

```{code-cell} ipython3
:tags: [remove-cell]
import msprime
import tskit

def sim_example1():
    ts = msprime.sim_ancestry(3, sequence_length=3e4,
            recombination_rate=1e-8, population_size=1000,
            record_full_arg=True, random_seed=21)
    ts.dump("data/haplotypes1.trees")

sim_example1()  # uncomment to recreate the tree seqs used in this notebook
```

```{code-cell} ipython3
ts = tskit.load("data/haplotypes1.trees")
ts.draw_svg(size=(400, 200), y_axis=True, time_scale='rank')
```

What we will do is to keep track of ancestrally inherited segments
as we move back through time, splitting and coalescing these as we go.
We will maintain the state as a dictionary mapping
each node ID to a list of segments.
What exactly does one of these lists of segments mean?
Well, a node represents an ancestral genome that was present
at a particular point in time.
An edge represents a sequence of ancestral genomes along which
a given segment was inherited.
(Anoter interpretation would be that an edge represents
a relationship: in other words, merely that a given node inherited
a given chunk of genome from a given ancestral node,
but let's go with the more concrete interpretation.)
So, to initialize we need to start off with a segment
that is the entire genome for each sample (since in this
tree sequence, all samples are at time=0 ago).

```{code-cell} ipython3
ancestry = {n.id : SegmentList() for n in ts.nodes()}
for j in ts.samples():
    ancestry[j].add_segment(0, ts.sequence_length)
```

Now, edges in the EdgeTable are sorted by parent time,
so if we iterate through the edges in order, we move back in time.
So, we can use this to see the state of the process at, say,
500 generations in the past:

```{code-cell} ipython3
for e in ts.edges():
    t = ts.node(e.parent).time
    if t > 300:
        break
    ca = ancestry[e.child]
    pa = ancestry[e.parent]
    seg = ca.remove_segment(e.left, e.right)
    for a, b in seg.segments:
        pa.add_segment(a, b)

for j in ancestry:
    print(j, ancestry[j].segments)
```

Here we see that the ancestral segments above nodes 2, 9 and 10
span the entire genome, while node 6 only spans part of the genome.
Since this tree sequence was simulated by msprime using
`record_full_arg=True`, interpretation is straightforward:
at 300 generations ago, there were three extant genomes
from which the samples inherited, and the inherited segments are
as listed here.
Note that this does not mean that "node 2 was laive 300 generations ago"
(clearly, as node 2 represetns an extant, sampled genome),
but rather that there are no other ancestral genomes recorded explicitly
in the tree sequence that lie on the path along
which node 2 has inherited it's genome.

We can use this to plot the state of the process as we move back through time.
Here are line segments depicting which bits of the genome are inherited
from which others at five times in the past.

```{code-cell} ipython3
import matplotlib.pyplot as plt
from matplotlib import collections as mc

times = [0, 50, 100, 200, 1000]
fig, axes = plt.subplots(len(times), figsize=(6, 8))

ancestry = {n.id : SegmentList() for n in ts.nodes()}
for j in ts.samples():
    ancestry[j].add_segment(0, ts.sequence_length)

k = 0
for e in ts.edges():
    t = ts.node(e.parent).time
    if t > times[k]:
        ax = axes[k]
        ax.set_title(f"t = {times[k]}")
        ax.set_xlim(0, ts.sequence_length)
        ax.set_ylim(0, ts.num_nodes)
        lc = mc.LineCollection(
            [(a, j), (b, j)] for j in ancestry
            for a, b in ancestry[j].segments
        )
        ax.add_collection(lc)
        k += 1
        if k == len(times):
            break
    ca = ancestry[e.child]
    pa = ancestry[e.parent]
    seg = ca.remove_segment(e.left, e.right)
    for a, b in seg.segments:
        pa.add_segment(a, b)

plt.tight_layout()
```

## Tracking ancestral segments

The code above gives an example of the basic algorithm that updates inherited segments.
However, for most purposes we want some additional structure.
For instance, let's say we'd like to know how many samples inherit from each segment.
To do this, we use a structure similar to the above,
but each segment now has a *label*, that is just the number of samples inheriting from it.
We'll demonstrate with a slightly more complex example tree sequence:

```{code-cell} ipython3
:tags: [remove-cell]
def sim_example2():
    ts = msprime.sim_ancestry(6, sequence_length=7e5,
            recombination_rate=1e-8, population_size=1000,
            random_seed=21)
    ts.dump("data/haplotypes2.trees")

sim_example2()  # uncomment to recreate the tree seqs used in this notebook
```

```{code-cell} ipython3
ts = tskit.load("data/haplotypes2.trees")
ts
```

Here is a data structure for a list of segments with labels:

```{code-cell} ipython3
class LabelSegmentList:
    
    def __init__(self, segments=None):
        if segments is None:
            segments = []
        self.segments = segments
    
    def __str__(self):
        return self.segments.__str__()
    
    def __repr__(self):
        return self.segments.__repr__()
    
    def __iter__(self):
        for abx in self.segments:
            yield abx
    
    def add_segment(self, left, right, label):
        self.segments.append([left, right, label])
        self.squash()
    
    def remove_segment(self, left, right):
        """
        Removes the intersection with [left, right).
        """
        removed = []
        new = []
        for a, b, x in self.segments:
            u, v = a, min(b, left)
            if u < v:
                new.append([u, v, x])
            u, v = max(a, left), min(b, right)
            if u < v:
                removed.append([u, v, x])
            u, v = max(a, right), b
            if u < v:
                new.append([u, v, x])
        self.segments = new
        return LabelSegmentList(removed)
    
    def iter_intersection(self, other):
        """
        Removes all portions of segments that intersect with any segment in
        other, iterating over theresulting removed segments, paired with the 
        label of the overlapping segment in other.
        """
        for left, right, label in other:
            intersection = self.remove_segment(left, right)
            for a, b, x in intersection.segments:
                yield [a, b, x], label
    
    def remove_intersection(self, other):
        """
        Simply removes the bits as in iter_intersection, but without returning them.
        """
        for _ in self.iter_intersection(other):
            pass
    
    def update_label(self, old, label):
        return old + label
    
    def squash(self):
        self.segments.sort()
        new = []
        i = 0
        while i < len(self.segments):
            a, b, x = self.segments[i]
            j = i + 1
            while (
                    j < len(self.segments)
                    and b == self.segments[j][0]
                    and x == self.segments[j][2]
                ):
                b = self.segments[j][1]
                j += 1
            new.append([a, b, x]) 
            i = j
        self.segments = new
    
    def merge(self, other):
        for (a, b, x), label in self.iter_intersection(other):
            new_x = self.update_label(x, label)
            self.segments.append([a, b, new_x])
        other.remove_intersection(self)
        for a, b, x in other:
            self.segments.append([a, b, x])
        self.squash()


```

For instance, now if we run this through to the end, 
we find all the roots (we can tell they are roots since their labels
are all "12", the number of samples), and which segments they are roots for:

```{code-cell} ipython3
ancestry = {n.id : LabelSegmentList() for n in ts.nodes()}
for j in ts.samples():
    ancestry[j].add_segment(0, ts.sequence_length, 1)

for e in ts.edges():
    t = ts.node(e.parent).time
    ca = ancestry[e.child]
    pa = ancestry[e.parent]
    seg = ca.remove_segment(e.left, e.right)
    pa.merge(seg)

for j in ancestry:
    print(j, ancestry[j])
```

Let's look at this in a different way:
here, each distinct segment at a few times
is plotted with *y*-coordinate equal to
the number of samples below that segment.
So, we look further back in time,
the segments break up into smaller pieces
and accumulate more samples as they start to coalesce.
(Segments are colored according to their node,
but not all are visible thanks to overplotting.)

```{code-cell} ipython3
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import colormaps
hsv = colormaps['hsv']

times = [0, 200, 400, 800, 1200, 2000, 4000]
fig, axes = plt.subplots(len(times), figsize=(6, 8), sharex=True)

ancestry = {n.id : LabelSegmentList() for n in ts.nodes()}
for j in ts.samples():
    ancestry[j].add_segment(0, ts.sequence_length, 1)

k = 0
for e in ts.edges():
    t = ts.node(e.parent).time
    if t > times[k]:
        ax = axes[k]
        ax.text(.01, .99, f"t = {times[k]}", ha='left', va='top', transform=ax.transAxes)
        ax.set_xlim(0, ts.sequence_length)
        ax.set_ylim(0, ts.num_samples + 1)
        # ax.set_ylabel("number of samples")
        lc = mc.LineCollection(
            [(a, x), (b, x)] for j in ancestry
            for a, b, x in ancestry[j].segments
        )
        lc.set_colors(list(hsv(j / ts.num_nodes) for j in ancestry
                      for _ in ancestry[j].segments))
        ax.add_collection(lc)
        k += 1
        if k == len(times):
            break
    ca = ancestry[e.child]
    pa = ancestry[e.parent]
    seg = ca.remove_segment(e.left, e.right)
    pa.merge(seg)

ax.set_xlabel("genomic position")
plt.tight_layout()
```

### Diversity calculation, moving back in time

We can use this data structure to compute mean tree distance
(i.e., twice the TMRCA) between the samples,
which is `ts.diversity(mode="branch")`.
To do this, we just need to keep track of how many distinct
pairs of lineages coalesced each time two segments merge.
To do this, we can just add two lines to the `merge( )` method
that returns this information.
We will then compute the mean TMRCA as
$$
    \pi = \frac{\sum_i t_i (b_i - a_i) n_i m_i}{ L n (n-1) / 2 },
$$
where the sum is over distinct events for which two distinct
lines of descent covering the segment $[a_i, b_i)$ coalesce;
the time this happened at is $t_i$,
the number of samples below each line of descent is $n_i$ and $m_i$,
the total sequence length is $L$, and the total number of samples is $n$.

```{code-cell} ipython3
class TmrcaLabelSegmentList(LabelSegmentList):
    
    def merge(self, other):
        pairs = 0
        for (a, b, x), label in self.iter_intersection(other):
            pairs += x * label * (b - a) # <--- new addition
            new_x = self.update_label(x, label)
            self.segments.append([a, b, new_x])
        other.remove_intersection(self)
        for a, b, x in other:
            self.segments.append([a, b, x])
        self.squash()
        return pairs # <--- and here

ancestry = {n.id : TmrcaLabelSegmentList() for n in ts.nodes()}
for j in ts.samples():
    ancestry[j].add_segment(0, ts.sequence_length, 1)

total = 0
for e in ts.edges():
    t = ts.node(e.parent).time
    ca = ancestry[e.child]
    pa = ancestry[e.parent]
    seg = ca.remove_segment(e.left, e.right)
    pairs = pa.merge(seg)
    total += pairs * t

total *= 2 / (ts.num_samples * (ts.num_samples - 1) * ts.sequence_length)

print(f"Branch-mode diversity: {ts.diversity(mode='branch')}")
print(f"   as calculated here: {2 * total}")
```

