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

Tree sequences can be output both by forwards-in-time genetic simulators
(such as SLiM) and forwards-in-time simulation libraries (such as
fwdpp/fwdpy11). If you are simply trying to obtain a tree sequence which is
the result of a forwards-in-time simulation, we recommend that you use one of
these tools, and consult the documentation and tutorials that go with them.

This tutorial *does not* attempt to show you how to use any existing forward
simulation tools. Instead, it shows the principles of using tskit to
build your own forwards-in-time tree sequence recording simulator from scratch.

:::{note}
Add details on building a forward simulator (see https://github.com/tskit-dev/tutorials/issues/14)
:::
