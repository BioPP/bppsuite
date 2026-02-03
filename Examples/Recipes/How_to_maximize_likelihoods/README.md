<!-- SPDX-FileCopyrightText: The Bio++ Development Group

     SPDX-License-Identifier: CECILL-2.1 -->


# Tutorial for first hints in maximum likelihood inference

Laurent Guéguen laurent.gueguen@univ-lyon1.fr


In this tutorial, we will see several examples of usage of programs
from bppsuite. Bppsuite programs depend on bio++ libraries, and all
are available
[there](https://github.com/BioPP/bpp-documentation/wiki).


We will give examples how to simulate sequences, perform maximum
likelihood estimates of models and branch lengths, perform ancestral
reconstruction and substitution mapping.

All the options are described in the bppsuite
[manual](https://pbil.univ-lyon1.fr/bpp-doc/bppsuite/bppsuite.html).



## Site process

The basic brick is the site process. A site process is the necessary
material to define the evolution of a site : a root, a tree, and a
model. In bio++, it is possible to mix several roots, trees and models
in many ways, to define simple to complex site process.


### Maximum Likelihood

In `hom_ML.bpp`, we define an homogeneous model on DNA evolution.
The results are stored in directory `hom_Res`, defined through
a macro. Make sure this directory exists before executing
`bppml`. Then the command is:

.```{bash}
bppml param=hom_ML.bpp
```

### Ancestral sequence inference

Given a process and data, it is possible to infer ancestral sequence
probabilities, using program \texttt{bppancestor}. Usually this is
done after the estimate of the best process given the data. In our
case, we use the parameter file `hom_Anc.bpp`.

```{bash}
bppancestor param=hom_Anc.bpp
```

### Non-homogeneous simple model

Now, we define several models of the tree, with specific assignations
to branches. For example, in file `nonhom_ML.bpp`.

In this case, two models are assigned to different branches. In the
case of trees read from newick files, the branches are numbered in the
reading order of the description string. In NHX format, branches are
explicitly numbered, which is much easier to handle.

### Aliases

We see in the declaration of second model

```
model2 = T92(theta=HKY85.theta_1)
```

that theta parameters of both models are linked, which means they are
constrained to share the same value (here equilibrium GC content).

Another way to alias parameters is through the line:

```
likelihood.alias = HKY85.kappa_1 -> T92.kappa_2
```

which means that kappa parameter model HKY85 is equal to kappa
parameter of model T92.

### Substitution rates

In this example, we define a distribution of substitution rates:

```
rate_distribution1 = Invariant(dist=Gamma(n=4))
```

which is here a Gamma distribution with 4 classes and a specific
probability of null rate for constant sites.

The posterior probabilities of the 5 classes of this
distribution are written in the info file `aln1.infos_1`. We can
see that for non-constant sites the posterior
probabilities of the first class (the invariant class) are quasi null.

### Computing on several data

#### Shared process

In `multidata_ML.bpp`, we compute the likelihood of the same
process as before on two alignments. 
In the result file `multiData_Res/aln1.params.txt`, the line

```
result=phylo1 + phylo2
```

means that the result (and optimized) log-likelihood is the sum of
log-likelihoods on both alignments. It means that the result is the
same as the log-likelihood of the concatenation of both alignments.

This result calculation is the default behaviour, but it is possible
to define its own calculation, such as: 

```
result=phylo1 + phylo2/2
```

or 

```
result=(phylo1 - phylo2) * (phylo2 - phylo1)
```

The optimization procedure is a minimization, so the formula must be
non positive.

#### Several processes

If we want to compute likelihood on several data with different trees,
we have to define several processes. See file `multiProc_ML.bpp`.

Both optimized trees are in the output directory. In this case, all
branch lengths are optimized independently. If we want to share branch
lengths between trees, we alias them. In this case, NHX format is more
convenient. See file `multiProcNHX_ML.bpp`, where the matching
branches of trees from files `aln1.nhx` and `aln2.nhx` are linked.

# Sequence processes

It is also possible to model the evolution of a sequence through
several processes, in a predetermined or probabilistic (mixture or
HMM) way. In the following examples, both processes differ on the
topology of the tree.

In `partition_ML.bpp`, an _a priori_ partition of the data
is proposed, defined in `process3`, which uses `process1`
and `process2`.

In `mixture_ML.bpp`, on each position the likelihood is the mean of
the likelihood of both processes, averaged through optimized
probabilities. To know which process best fits each site, in
`mixture_Res/infos`, the site log-likelihoods and posterior
probabilities are displayed.

In `autocorr_ML.bpp`, the modeling is also a probabilistic mixture of
both models, but with a between sites auto-correlation for each
process.


\end{document}
