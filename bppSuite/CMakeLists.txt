# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 22/08/2009

add_executable(bppml bppML.cpp)
target_link_libraries(bppml ${LIBS})

add_executable(bppseqgen bppSeqGen.cpp)
target_link_libraries(bppseqgen ${LIBS})

add_executable(bppdist bppDist.cpp)
target_link_libraries(bppdist ${LIBS})

add_executable(bpppars bppPars.cpp)
target_link_libraries(bpppars ${LIBS})

add_executable(bppseqman bppSeqMan.cpp)
target_link_libraries(bppseqman ${LIBS})

add_executable(bppconsense bppConsense.cpp)
target_link_libraries(bppconsense ${LIBS})

add_executable(bppancestor bppAncestor.cpp)
target_link_libraries(bppancestor ${LIBS})

add_executable(bppreroot bppReRoot.cpp)
target_link_libraries(bppreroot ${LIBS})

add_executable(bppphysamp bppPhyloSampler.cpp)
target_link_libraries(bppphysamp ${LIBS})

add_executable(bpptreedraw bppTreeDraw.cpp)
target_link_libraries(bpptreedraw ${LIBS})

# Install libs
install(TARGETS
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
