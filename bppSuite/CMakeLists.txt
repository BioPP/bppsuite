# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 22/08/2009

ADD_EXECUTABLE(bppml bppML.cpp)
TARGET_LINK_LIBRARIES(bppml ${LIBS})

ADD_EXECUTABLE(bppseqgen bppSeqGen.cpp)
TARGET_LINK_LIBRARIES(bppseqgen ${LIBS})

ADD_EXECUTABLE(bppdist bppDist.cpp)
TARGET_LINK_LIBRARIES(bppdist ${LIBS})

ADD_EXECUTABLE(bpppars bppPars.cpp)
TARGET_LINK_LIBRARIES(bpppars ${LIBS})

ADD_EXECUTABLE(bppseqman bppSeqMan.cpp)
TARGET_LINK_LIBRARIES(bppseqman ${LIBS})

ADD_EXECUTABLE(bppconsense bppConsense.cpp)
TARGET_LINK_LIBRARIES(bppconsense ${LIBS})

ADD_EXECUTABLE(bppancestor bppAncestor.cpp)
TARGET_LINK_LIBRARIES(bppancestor ${LIBS})

ADD_EXECUTABLE(bppreroot bppReRoot.cpp)
TARGET_LINK_LIBRARIES(bppreroot ${LIBS})

ADD_EXECUTABLE(bppphysamp bppPhyloSampler.cpp)
TARGET_LINK_LIBRARIES(bppphysamp ${LIBS})

ADD_EXECUTABLE(bpptreedraw bppTreeDraw.cpp)
TARGET_LINK_LIBRARIES(bpptreedraw ${LIBS})

# Install libs
INSTALL(TARGETS
  bppml
  bppseqgen
  bppdist
  bpppars
  bppseqman
  bppconsense
  bppancestor
  bppreroot
  bppphysamp
  bpptreedraw
  DESTINATION bin)
