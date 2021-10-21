### Introduction

This is an implementation of the MBC and MBC* algorithm written on the paper:

*Lyu, Bingqing, et al. "Maximum biclique search at billion scale." Proceedings of the VLDB Endowment (2020).*

It is written in C with a Python interface.

### Definition

Given a bipartite graph $$G = (U, V, E)$$, a **biclique** $$C$$ is a complete bipartite subgraph of
$$G$$, i.e., for each pair of $$u \in U(C)$$ and $$v \in V(C)$$, we have $$(u, v) \in E(C)$$.

In this implementation, we aim to find a biclique $$C^\*$$ in $$G$$ with the maximum number of **edge**.

If you are looking for **enumerating** maximum biclique, I suggest you refer to:

*Alexe, Gabriela, et al. "Consensus algorithms for the generation of all maximal bicliques." Discrete Applied Mathematics 145.1 (2004): 11-21.*

*Lu, Yuping, Charles A. Phillips, and Michael A. Langston. "Biclique: an R package for maximal biclique enumeration in bipartite graphs." BMC research notes 13.1 (2020): 1-5.*

and the implementation [here](https://cran.r-project.org/web/packages/biclique/index.html).

If you are looking for finding the biclique with the maximum number of **vertex**, I suggest you refer to the [NetworkX](https://networkx.org/) project and this [post](https://cs.stackexchange.com/questions/131081/polynomial-time-algorithm-to-solve-the-maximum-vertex-bipartite-subgraph-problem).

### About the implementation

When I implemented the algorithm, I had an application in my mind. As a result, it probably does not suitable for your own use. For example, D is passed as a dense matrix and the bipartite graph is represented as two arrays of bitset internally to speed up some steps. 

I also used `_GNU_SOURCE` to use some `stdlib.h` method. It may limit the portability of the code.

Good luck.

### Installation

#### Prerequisites
You need [setuptools](https://pypi.org/project/setuptools/), https://numpy.org/[numpy](numpy) and standard build toolchain. If you are running on a Ubuntu, the following command should suffice:

```bash
$ apt install python3 python3-numpy python3-setuptools build-essential libpython3-dev
```

#### Compile and Installation

```bash
$ python setup.py build
$ sudo python setup.py install # or python setup.py install --user
```

### Documentstion

```python
import pymbc
help(pymbc)
```

### Examples

Please refer to the `examples/` subdirectory.

### License

3-Clause BSD License

### Contact

WONG Hang wonghang (at) gmail (dot) com
WONG is my surname.

