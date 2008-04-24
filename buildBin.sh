#! /bin/sh

strip bppSuite/bppdist
strip bppSuite/bpppars
strip bppSuite/bppml
strip bppSuite/bppseqgen
strip bppSuite/bppconsense
strip bppSuite/bppseqman
strip bppSuite/bppphysamp
strip bppSuite/bppreroot
tar cvzf bppSuite-linux32-bin-static-0.2.3.tar.gz bppSuite/bppdist bppSuite/bpppars bppSuite/bppml bppSuite/bppseqgen bppSuite/bppconsense bppSuite/bppseqman bppSuite/bppphysamp bppSuite/bppreroot

