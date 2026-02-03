# SPDX-FileCopyrightText: The Bio++ Development Group
#
# SPDX-License-Identifier: CECILL-2.1

#! /bin/sh
arch=`uname -m`
version=3.0.0-1

strip bppSuite/bppdist
strip bppSuite/bpppars
strip bppSuite/bppml
strip bppSuite/bppseqgen
strip bppSuite/bppconsense
strip bppSuite/bppseqman
strip bppSuite/bppreroot
strip bppSuite/bppancestor
strip bppSuite/bpptreedraw
strip bppSuite/bppalnscore
strip bppSuite/bppmixedlikelihoods
strip bppSuite/bpppopstats
strip bppSuite/bppbranchlik
tar cvzf bppsuite-${arch}-bin-static-${version}.tar.gz bppSuite/bppdist bppSuite/bpppars bppSuite/bppml bppSuite/bppseqgen bppSuite/bppconsense bppSuite/bppseqman bppSuite/bppreroot bppSuite/bppancestor bppSuite/bpptreedraw bppSuite/bppalnscore bppSuite/bppmixedlikelihoods bppSuite/bpppopstats bppSuite/bppbranchlik

