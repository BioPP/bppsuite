#! /bin/sh

strip bppSuite/bppdist
strip bppSuite/bpppars
strip bppSuite/bppml
strip bppSuite/bppseqgen
strip bppSuite/bppconsense
tar cvzf bppSuite-linux-bin-static-0.1.0.tar.gz bppSuite/bppdist bppSuite/bpppars bppSuite/bppml bppSuite/bppseqgen bppSuite/bppconsense

