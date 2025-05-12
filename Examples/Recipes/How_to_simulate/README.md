<!-- SPDX-FileCopyrightText: The Bio++ Development Group

     SPDX-License-Identifier: CECILL-2.1 -->

# Tutorial for bppsuite

Laurent Guéguen <laurent.gueguen@univ-lyon1.fr>

In this tutorial, we will see how simulation of sequence evolution can
be performed using the Bio++ libraries and the `bppseqgen` program, on
which all information is available [here](https://github.com/BioPP/bpp-documentation/wiki).

All the options are described in the [bppsuite manual](https://pbil.univ-lyon1.fr/bpp-doc/bppsuite/bppsuite.html).

We will specifically explore the simulation of protein alignments, but
these approaches are also available for nucleotide and codon
alignments.

A simulation follows a subtitution process, which is defined by
a (set of) substitution models, a tree, and optionally, a root
distribution and a substitution rate distribution.

- In all simulations, a tree is needed to guide the evolution
  process. We first consider the tree in file `../../Data/LSUrooted.dnd`.
- Many amino acid models are available in Bio++, listed [here](https://pbil.univ-lyon1.fr/bpp-doc/bppsuite/bppsuite.html#Protein).
  
In the next sections, we will see increasingly complex ways to perform simulations
with such models, allowing to extract more and more information from the data.


## Homogeneous in time and space

First, we describe examples where the model defines the same process on all
branches and sites.

### Simple homogeneous model

In `simple_hom.bpp`, we consider the simplest evolutionary scenario:
homogeneous and stationnary, with a simple substitution model
 [WAG01](https://pbil.univ-lyon1.fr/bpp-doc/bpp-phyl/html/classbpp_1_1WAG01.html#details),
following the tree in file `LSU.dnd`. Both are gathered into
process `process1`. `simul1` sets the simulation, using
this process, of a 300 amino acids protein alignment.

```{bash}
bppseqgen param=simple_hom.bpp
```

In `simul1`, the argument `output.internal.sequences = true` means that the output file also contains the simulated ancestral sequences, otherwise
only sequences at leaf nodes (extant species) are included. The sequences are labeled according to the
node ids.

### Homogeneous mixture model

A biology-realistic model consists of a mixture of models, where sites can evolve according to different profiles. In the next examples, we use the
[LGL08_CAT](https://pbil.univ-lyon1.fr/bpp-doc/bpp-phyl/html/classbpp_1_1LGL08__CAT.html#details) family of models
that uses a given number of amino acids profiles as found in databases.

In addition, we introduce a variable (NBCAT) that can be set from the
command line:

```{bash}
bppseqgen param=mixed_hom.bpp NBCAT=40
```

In this file, two simulations are performed. 

There are many ways to organize the mixing of models in a tree,
depending on how the choice of a model on a branch depends on the choice
on the upper branch. In file `mixed_hom.bpp`, the line

```
scenario1=split(model=1)
```

defines a simple modle *organizaiton scenario*: a
site will follow the same model all along the tree, which with CAT
means the same profile. `process1` follows this scenario, for
`simul1`.

Without this option, at each node the process can switch freely
between the submodels of the mixture. `process2` defines such a
simulation, and is used in simulation `simul2`.

Much more elaborate scenarios can be defined with the definition of
`paths` (see below?)


### Substitution rates

In the previous examples, all sites evolved at the same rate, whereas
we could wish an heterogeneity in these rates. This can be done in two
ways, either in a global or site-specific (see section
[Non-homogeneous in space](#sec-non-homog-space) manner.

In this example, we define a distribution of substitution rates:

```
rate_distribution1 = Invariant(p=0.15, dist=Gamma(n=4))
```

which is, here, a discrete Gamma distribution with 4 classes and a
probability 0.15 of constant sites (null rate).


### Root frequencies

In the previous descriptions, the modeling is stationary, which means
that the states distribution at the root is the stationary
distribution of the model.

It is possible to set specific root frequencies (or states, see below)
to define non-stationary processes. In file `mixed_IG_hom.bpp`, a
`root_freq` is defined, and used in the definition of
`process2` as the root frequencies, which will be used for
`simul2`.

Both simulations can be run with the command:

```
bppseqgen param=mixed_IG_hom.bpp NBCAT=40
```


## <a name="sec-non-homog-space"></a>Non-homogeneous in space

In the previous models, the process is defined similarly for all
sites (even though with mixture models a site-specific choice of
submodels is achieved). It is also possible to define site-specific models, for the
root states, the substitution rates, and the substitution models.

This is independent of modeling time non-homogeneity, so both can be
used together. 


### Site-specific root states

States at the root of the tree can be set or sampled. In the tsv file
`infos_rates.tsv`, column `States` assigns states to sites
at the root node. In file `sites_rate_state_hom.bpp`, this information is used
for `simul1`. The column of states can be named differently, as
long as it is consistent with `input.infos.states`.

Another way to define states at the root is to use a squence alignment, specified by a
given file path. Still in `sites_rate_state_hom.bpp`,
`input.data1` is an alignment read from file `LSU.phy`,
from which a random sequence is used as a root in `simul2`.

If the simululation uses a process with a defined `root_freq`, the states
defined with `input.infos.states` are prioritary (EDIT[Julien]: I do not understand this part...).


### Site-specific substitution rates

Rates can also set up in a site-specific way, in a similar tsv file.
For example in column `Rates` of file `infos_rates.tsv`,
the middle sites are slower than the extreme sites, as used in `simul3`
of `sites_rate_state_hom.bpp`.

If the simulation uses a process with a defined `rate`, the rates
defined with `input.infos.rates` are prioritary.


The rates can be defined without the root states (in which case those
states are sampled from a distribution). In `simul4`, they are
taken from a sequence in `input.data1`. In this case, only the
first 981 sites are considered (to fit with the number of rates in
`infos_rates.tsv`). 

```{bash}
bppseqgen param=sites_state_rate_hom.bpp NBCAT=40
```

still with the site-specific rates defined in `infos_rates.tsv`.

### Site-specific processes

Beyond site-specific rates and root states, it is also possible to set
up site-specific processes. In the `partition.bpp` file, two simple
processes are defined, following the same stationary WAG01 model but
different trees. Then they are assigned to different sets of sites, in
a new process, which is used for simulation.

```{bash}
bppseqgen param=partition.bpp
```

Instead of presetting boundaries in the partition, processes can be put
into a markovian chain along the sequence, or more simply with
autocorrelation probabilities of successive presence, as in
`autocorr.bpp` for the same processes as before.

## Non-homogeneous in time

In addition to the non-homogeneity in space, it is possible, totally
independently, to model time heterogenity, which means that
different models can be assigned to different branches of the tree.

### Assigned non-homogeneity
  
Still in file `nonstat_WAG_nonhom.bpp`, `process1` is a
non-homogeneous process made of two WAG01 with
different equilibrium frequencies, assgined to two sets of branches,
and specific root frequencies.

`simul1` performs simple simulations with this model, whereas
`simul2` also uses site-specific rates

### Non-homogeneous mixtures

This can be done also with mixture models, and all possibilities of
usage of submodels can be set within scenarios. In
`nonstat_CAT_nonhom.bpp`, two CAT models are used, and with the
description:

```
scenario1=split(model=(1,2))
```

all paths where an identical submodel of each CAT model is used on all
branches where this model is set. For example, with command:

```{bash}
bppseqgen param=nonstat_CAT_nonhom.bpp NBCAT=10
```

all 100 possible paths are with 2 submodels of CAT 10 are possible.

## Posterior simulations

In all previous examples, simulations were done using given models, root
frequencies, and branch lengths. But in the prospect of simulations as
realistic as possible, it is difficult to choose without information about
real data. In Bio++, all those parameters can be estimated through
maximum likelihood inference (via `bppml`), and given as inputs
for simulations. But because of the numerous hypotheses set in the
modeling (such as site-independence, too simple models, etc), the
resulting simulations can still be very unrealistic.

Another way to include data information in the simulation is to use
posterior transition probabilities. On a branch and a site, if the
used model is $M$ and the data $D$, for two states $x$ and $y$,
instead of using $P(y|x,M)$ as transition probabilities, we use
$P(y|x,M,D)$.

Hence site and branch specific information is taken into account, even
if it is not explicitly modeled.

In Bio++, the `phylo` object links a `process` with
`data`, providing a likelihood. To simulate alignments
*a posteriori*, one just need to replacethe reference
to `process` with a reference to `phylo` in the `simul` description.

In a complete posterior simulation, the sequences simulated at the
leaves will be exactly those of the given alignment, which is not very
interesting. So it is possible to release this to use directly the
models (and the prior probabilities) on some branches. In
`simul`, the `nullnodes` argument accepts a list of branch
ids, where the data will not be taken into account.

In the `simple_data.bpp` file, `nullnodes` is set to
shortcut `Leaves`, which means that the simulation on all
external branches is following the model.

```{bash}
bppseqgen param=simple_data.bpp
```

Beware that such simulations takes more time than previouly, since all
the conditional likelihoods have to be computed.

