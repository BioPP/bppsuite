<!-- SPDX-FileCopyrightText: The Bio++ Development Group

     SPDX-License-Identifier: CECILL-2.1 -->

# Tutorial for infering GC3*

Laurent Guéguen <laurent.gueguen@univ-lyon1.fr>

In this tutorial, we will see how Bio++ libraries can be used to
explore the heterogeneity of equilibrium GC3 (aka GC3*) in codon
alignments, among branches. All information about Bio++ is available
[here](https://github.com/BioPP/bpp-documentation/wiki).

In addition, all the options are described in the [bppsuite
manual](https://pbil.univ-lyon1.fr/bpp-doc/bppsuite/bppsuite.html)
and [testnh manual](https://pbil.univ-lyon1.fr/bpp-doc/testnh/testnh.html).


We will proceed through homogeneous and branch-specific modeling.

## YN98 

### Homogeneous modeling

First, the model defines the same process on all branches and sites.

#### Maximum likelihood

In `hom_YN98_ML.bpp`, we consider a homogeneous process, using a simple codon substitution model
[YN98](https://pbil.univ-lyon1.fr/bpp-doc/bpp-phyl/html/classbpp_1_1YN98.html#details),
following the tree in file `spec_tree.dnd`. The process is denoted
`process1`. 

Stationary distribution of any codon is defined through the product of
the frequencies of its position specific nucleotides (F3X4). The
frequencies are normalized after removing stop codons. Other codon
frequencies are available, see
[there](https://pbil.univ-lyon1.fr/bpp-doc/bppsuite/bppsuite.html#Codon-alphabets).

The results are stored in directory `hom_Res`, defined through a
macro. Make sure this directory exists before executing
\texttt{bppml}. Then the command is:

```{bash}
  bppml param=hom_YN98_ML.bpp
```

A backup file `hom_Res/params.bup` lists all parameters values as they
are optimized, and is used as restarting point in case of issue. 

In `hom_Res/params.txt`, the most likely modelk is displayed, using
the same syntax as in `.bpp` file. The most likely tree is in
`hom_Res/tree_ML.dnd_1`.


#### Substitution Mapping 

Given a process and data, the program `mapnh` (from
[Testnh](https://github.com/BioPP/testnh) computes the expected number
of ancestral events of specific types on sites and/or branches.

Any type of events can be counted, and here we are interested counting
the substitutions between AT3 and GC3 on all branches.

```{bash}
  mapnh param=map_GC3.bpp
```

The output is a table with branch numbers in columns and substitution
types in rows. Without the option `format=tsv` in `output.counts`, the
output is in a newick format.


The previous mapping computes raw counts, from which the trend towards
GC3 can be computed as $\frac{\text{GC3}}{\text{AT3 } + \text{ GC3}}$.
But this ratio does not account from the ancestral sequence
composition, which can bias this ratio. To fix this, we compute the
expected number of same events that could have been performed with a
model where the equilibrium GC3 is balanced, 0.5. We call these values
as **opportunities**.

In the previous model, GC3 is modeled through parameter
`YN98.3_Full.theta_1` in F3X4 (details 
[here](https://pbil.univ-lyon1.fr/bpp-doc/bpp-phyl/html/classbpp_1_1CodonFromIndependentFrequencySet.html#details))


In `map_GC3_norm.bpp`, the option `nullProcessParams` allows this
computation, setting `YN98.3_Full.theta_*=0.5` (with * as a wildcard)
for the GC3-unbiased model.

```{bash}
  mapnh param=map_GC3_norm.bpp
```

Two files are produced, for counts and for opportunities. The
normalized counts, per unit time, are then the ratio counts over
opportunities, from which the GC3 trend is computed as before. Without
option ``splitNorm=true'', the normalized ratios are directly output.


### Maximum likelihood estimates of branch GC3

It is possible to explore GC3* modeling heterogeneous evolution in the
tree A way to do this is to optimize one model per branch. In this
case, there would be 12 parameters per branch, in addition to the root
frequencies, which is too many. So we consider some assumptions:

1. kappa and omega, equilibrium frequencies excepted GC3 are shared
   among branches; For this we will use the `shared_parameters` option
   in the description of the process;
2. branch lengths are not re-estimated from homogeneous modeling; For
   this we use option `ignore_parameters` in the optimization options
   and we use the tree formerly estimated in the homogeneous modeling;
3. The optimized homogeneous process is a starting point for the
   optimization process.


```{bash}
bppml param=opb_ML.bpp
```
The results are stored in directory `opb_Res`, defined through
a macro. Make sure this directory exists before.


The optimized modeling is output in file ``opb_Res/params.txt'': all
models are similar excepted ``3_Full.theta''. Those values can also be
found in ``opb_Res/params.bup.def''.


This parameterization gives an example of options of optimization that
can be used with BppML, but they can be removed (at the cost of much
much longer optimization).


