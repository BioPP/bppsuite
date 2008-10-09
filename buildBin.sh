#! /bin/sh

strip bppSuite/bppdist
strip bppSuite/bpppars
strip bppSuite/bppml
strip bppSuite/bppseqgen
strip bppSuite/bppconsense
strip bppSuite/bppseqman
strip bppSuite/bppphysamp
strip bppSuite/bppreroot
strip bppSuite/bppancestor
tar cvzf bppSuite-linux64-bin-static-0.3.0.tar.gz bppSuite/bppdist bppSuite/bpppars bppSuite/bppml bppSuite/bppseqgen bppSuite/bppconsense bppSuite/bppseqman bppSuite/bppphysamp bppSuite/bppreroot bppSuite/bppancestor

