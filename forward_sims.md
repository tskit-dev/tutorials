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

(sec_tskit_forward_simulations)=

# _Building a forward simulator_

% remove underscores in title when tutorial is complete or near-complete

This tutorial shows how tskit can be used to
build your own forwards-in-time tree sequence recording simulator from scratch.

:::{note}
If you are simply trying to obtain a tree sequence which is
the result of a forwards-in-time simulation, this can be done by using one of the
highly capable forwards-in-time genetic simulators that already exist, such as
[SLiM](https://messerlab.org/slim/) or [fwdpy11](https://github.com/molpopgen/fwdpy11).
Documentation and tutorials for these tools exist on their respective websites. This
tutorial focusses instead on illustrating the general principles that lie behind such
simulators.
:::

:::{todo}
Add details on building a forward simulator (see issue
[#14](https://github.com/tskit-dev/tutorials/issues/14))
:::
