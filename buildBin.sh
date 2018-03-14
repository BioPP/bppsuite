#! /bin/sh
arch=x86_64 #i686
version=2.3.2

strip bppSuite/bppdist
strip bppSuite/bpppars
strip bppSuite/bppml
strip bppSuite/bppseqgen
strip bppSuite/bppconsense
strip bppSuite/bppseqman
strip bppSuite/bppphysamp
strip bppSuite/bppreroot
strip bppSuite/bppancestor
strip bppSuite/bpptreedraw
strip bppSuite/bppalnscore
strip bppSuite/bppmixedlikelihoods
strip bppSuite/bpppopstats
tar cvzf bppsuite-${arch}-bin-static-${version}.tar.gz bppSuite/bppdist bppSuite/bpppars bppSuite/bppml bppSuite/bppseqgen bppSuite/bppconsense bppSuite/bppseqman bppSuite/bppphysamp bppSuite/bppreroot bppSuite/bppancestor bppSuite/bpptreedraw bppSuite/bppalnscore bppSuite/bppmixedlikelihoods bppSuite/bpppopstats

