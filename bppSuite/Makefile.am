#Copyright or © or Copr. CNRS, (November 16, 2004-2008)
#
#This software is a computer program whose purpose is to provide tools
#for sequence analysis, phylogenetics and molecular evolution.
#
#This software is governed by the CeCILL  license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL
#license as circulated by CEA, CNRS and INRIA at the following URL
#"http://www.cecill.info". 
#
#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 
#
#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 
#
#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL license and that you accept its terms.

ML_H_FILES = 

ML_CPP_FILES = \
bppML.cpp

DIST_H_FILES =

DIST_CPP_FILES = \
bppDist.cpp

PARS_H_FILES = 

PARS_CPP_FILES = \
bppPars.cpp

SEQGEN_H_FILES = 

SEQGEN_CPP_FILES = \
bppSeqGen.cpp

CONSENSE_H_FILES = 

CONSENSE_CPP_FILES = \
bppConsense.cpp

SEQMAN_H_FILES = 

SEQMAN_CPP_FILES = \
bppSeqMan.cpp

PHYSAMP_H_FILES = 

PHYSAMP_CPP_FILES = \
bppPhyloSampler.cpp

REROOT_H_FILES = 

REROOT_CPP_FILES = \
bppReRoot.cpp

AUTOMAKE_OPTIONS = gnu

bin_PROGRAMS = bppml bpppars bppdist bppseqgen bppconsense bppseqman bppphysamp bppreroot
if STATIC
bppml_LDFLAGS = -static
bpppars_LDFLAGS = -static
bppdist_LDFLAGS = -static
bppseqgen_LDFLAGS = -static
bppconsense_LDFLAGS = -static
bppseqman_LDFLAGS = -static
bppphysamp_LDFLAGS = -static
bppreroot_LDFLAGS = -static
else
bppml_LDFLAGS = 
bpppars_LDFLAGS = 
bppdist_LDFLAGS =
bppseqgen_LDFLAGS =
bppconsense_LDFLAGS =
bppseqman_LDFLAGS =
bppphysamp_LDFLAGS =
bppreroot_LDFLAGS =
endif
bppml_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq -lbpp-phyl
bppml_SOURCES = $(ML_H_FILES) $(ML_CPP_FILES)
bpppars_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq -lbpp-phyl
bpppars_SOURCES = $(PARS_H_FILES) $(PARS_CPP_FILES)
bppdist_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq -lbpp-phyl
bppdist_SOURCES = $(DIST_H_FILES) $(DIST_CPP_FILES)
bppseqgen_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq -lbpp-phyl
bppseqgen_SOURCES = $(SEQGEN_H_FILES) $(SEQGEN_CPP_FILES)
bppconsense_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq -lbpp-phyl
bppconsense_SOURCES = $(CONSENSE_H_FILES) $(CONSENSE_CPP_FILES)
bppseqman_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq
bppseqman_SOURCES = $(SEQMAN_H_FILES) $(SEQMAN_CPP_FILES)
bppphysamp_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq -lbpp-phyl
bppphysamp_SOURCES = $(PHYSAMP_H_FILES) $(PHYSAMP_CPP_FILES)
bppreroot_LDADD = -lbpp-utils -lbpp-numcalc -lbpp-seq -lbpp-phyl
bppreroot_SOURCES = $(REROOT_H_FILES) $(REROOT_CPP_FILES)

