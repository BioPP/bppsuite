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

AUTOMAKE_OPTIONS = gnu

bin_PROGRAMS = bppml bpppars bppdist bppseqgen bppconsense
if STATIC
bppml_LDFLAGS = -static
bpppars_LDFLAGS = -static
bppdist_LDFLAGS = -static
bppseqgen_LDFLAGS = -static
bppconsense_LDFLAGS = -static
else
bppml_LDFLAGS = 
bpppars_LDFLAGS = 
bppdist_LDFLAGS =
bppseqgen_LDFLAGS =
bppconsense_LDFLAGS =
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

