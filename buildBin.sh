#! /bin/sh
arch=linux64
version=0.4.0

strip bppSuite/bppdist
strip bppSuite/bpppars
strip bppSuite/bppml
strip bppSuite/bppseqgen
strip bppSuite/bppconsense
strip bppSuite/bppseqman
strip bppSuite/bppphysamp
strip bppSuite/bppreroot
strip bppSuite/bppancestor
tar cvzf bppSuite-${arch}-bin-static-${version}.tar.gz bppSuite/bppdist bppSuite/bpppars bppSuite/bppml bppSuite/bppseqgen bppSuite/bppconsense bppSuite/bppseqman bppSuite/bppphysamp bppSuite/bppreroot bppSuite/bppancestor

