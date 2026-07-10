<!-- SPDX-FileCopyrightText: The Bio++ Development Group

     SPDX-License-Identifier: CECILL-2.1 -->

# Tutorial for exploring site positive selection

Laurent Guéguen <laurent.gueguen@univ-lyon1.fr>

In this tutorial, we will see how Bio++ libraries can be used to
explore site-specific positive selection. All information about Bio++
is available [here](https://github.com/BioPP/bpp-documentation/wiki).
In addition, the options are described in the [bppsuite
manual](https://pbil.univ-lyon1.fr/bpp-doc/bppsuite/bppsuite.html).


The point of this tutorial is to present how successive models and
methods can be used; the hypothesis for the optimization process can
be questionned.

Our study will be upon primate sequences of APOBEC3G gene.

## Homogeneous in time and space

First, we describe examples where the model defines the same process
on all branches and sites.

### Simple homogeneous model

In `hom_M0.bpp`, we consider a homogeneous process, using a simple codon substitution model
[YN98](https://pbil.univ-lyon1.fr/bpp-doc/bpp-phyl/html/classbpp_1_1YN98.html#details),
following the tree in file `APOBEC3G_primates.dnd`. The process is denoted
`process1`. 

We used a Gamma+Invariant distribution of susbtitution rates: 

```
rate_distribution1 = Invariant(p=0.15, dist=Gamma(n=4))
```

which is, here, a discrete Gamma distribution with 4 classes and a
probability 0.15 of constant sites (null rate).


The results are stored in directory `hom_M0`, defined through
a macro. Make sure this directory exists before executing
bppml. Then the command is:

```{bash}
  bppml param=hom_M0.bpp
```

A backup file `hom_M0/APOBEC3G_primates_params.bup` lists all parameters values as they
are optimized, and is used as restarting point in case of issue. 

In `hom_M0/APOBEC3G_primates_params.bpp`, the most likely output is displayed, using
the same syntax as in `.bpp` file.


### Extend to homogeneous Mixture models

A usual way to explore site heterogeneity is to model the alignment
through mixture models. In Bio++, many kinds of mixtures are
available, but in the case of YN98 for the detection of selection, the
usual codeml models M1, M2, ... are available, and denoted
YNGP\_M....

Those models can be infered from scratch, but we show how to use
previous optimization to avoid useless (and messy) optimization.

For M1 modeling, the results are stored in directory `hom_M1`. Also,
make sure this directory exists before executing bppml. Likewise for
M2 model. 

1. For M1, in file `hom_M1.bpp`, `param` reads the modeling obtained in
   M0 optimization in file `hom_M0/APOBEC3G_primates_params.bpp`,
   and optimized tree is used as input tree.
   
   From this, all parameters infered in the first process are frozen
   with wildcard `*_1` in `ignore_parameters` configuration line.
   
   The M1 modeling is defined through `model2`, which is then called
   in `process1` and then `phylo1` for likelihood calculation.
   
2. branch lengths are not re-estimated; Argument `BrLen` in 
   `ignore_parameters` would do this, but there are already considered
   with `*_1` wildcard;
3. The optimized M0 homogeneous process is a starting point for the
   optimization process.


Between M1 and M2 modelings, the numbers of processes are different,
to allow joint usage of those processes.

```{bash}
bppml param=hom_M1.bpp

bppml param=hom_M2.bpp
```

### Homogeneous site analysis

In files `APOBEC3G_primates_infos.txt_1`, we can get per site log-likelihoods:

A quick glimpse in R:

```{R}
infosM1 = read.table("hom_M1/APOBEC3G_primates.ml_h.infos_1",h=T)
infosM2 = read.table("hom_M2/APOBEC3G_primates.ml_h.infos_1",h=T)

plot(infosM2$lnL-infosM1$lnL,t="l")
abline(h=0)

```

On mixture models, more information can be obtained with
`bppmixedlikelihoods` program, that computes log-likelihoods and
posterior probabilities of the sub-models of the mixture:

```{bash}
bppmixedlikelihoods param=mixedlikelihoods.bpp,hom_M2/APOBEC3G_primates.ml_h.params.bpp OUTPUTDIR=hom_M2
```

And if we want to plot the mean posterior omega with difference in
log-likelihood:

```{R}
llM2 = read.table("hom_M2/APOBEC3G_primates_llik.tsv",h=T)

plot(infosM2$lnL-infosM1$lnL,llM2$mean)

```

Also, branch-site information is available with program
`bppbranchlik`. Phylo-likelihood is computed with a given process, and
an alternative process is used iteratively on each branch to compute
the log-likelihood. 

In this example, the M1 model is used as background model (used for
the computation of `phylo1`), and M2 (in process 3, with
`alt_process=3`) is the foreground model.

```{bash}
bppbranchlik param=branchlik.bpp,hom_M1/APOBEC3G_primates.ml_h.params.bpp,hom_M2/APOBEC3G_primates.ml_h.params.bpp
```

And to plot:

```{R}
llM1M2 = read.table("APOBEC3G_primates_M1M2.tsv",h=T)
plot(llM1M2$Llik)

```
The resulting likelihoods can then be compared to the likelihood of
the M1 basic process to get 'nearly' significance (since this is not
the optimized model per branch).

THe same analysis can be done on branch-site basis.


## Clade analysis

In the next example, we discriminate a clade from the rest of the
tree, with M1 modeling against M2 modeling in this specific clade,
with a M0 background. Build directories `clade_M0M1` and `clase_M0M2`
before:


```{bash}
bppml param=clade_M0M1.bpp

bppml param=clade_M0M2.bpp
```

And then for site specific likelihoods that focus on this clade,
likewise:


```{R}
infosM1 = read.table("clade_M0M1/APOBEC3G_primates_infos.txt_1",h=T)
infosM2 = read.table("clade_M0M2/APOBEC3G_primates_infos.txt_1",h=T)

plot(infosM2$lnL-infosM1$lnL,t="l")
abline(h=0)

```
