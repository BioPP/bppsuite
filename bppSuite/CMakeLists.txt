# CMake script for Bio++ Program Suite
# Authors:
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 22/08/2009

# Executables of bppsuite.
# Generation of targets from file name is not automated in case of executables not following the pattern.

add_library(bppTools bppTools.cpp)

add_executable (bppml bppML.cpp)
add_executable (bppseqgen bppSeqGen.cpp)
add_executable (bppdist bppDist.cpp)
add_executable (bpppars bppPars.cpp)
add_executable (bppseqman bppSeqMan.cpp)
add_executable (bppconsense bppConsense.cpp)
add_executable (bppancestor bppAncestor.cpp)
add_executable (bppmixedlikelihoods bppMixedLikelihoods.cpp)
add_executable (bppreroot bppReRoot.cpp)
add_executable (bpptreedraw bppTreeDraw.cpp)
add_executable (bppalnscore bppAlnScore.cpp)
add_executable (bpppopstats bppPopStats.cpp)

set (bppsuite-targets
  bppml
  bppseqgen
  bppdist
  bpppars
  bppseqman
  bppconsense
  bppancestor
  bppmixedlikelihoods
  bppreroot
  bpptreedraw
  bppalnscore
  bpppopstats
  )

if (BUILD_STATIC)
  target_link_libraries (bppTools ${BPP_LIBS_STATIC})
  set_target_properties (bppTools PROPERTIES LINK_SEARCH_END_STATIC TRUE)
else (BUILD_STATIC)
  target_link_libraries (bppTools ${BPP_LIBS_SHARED})
  set_target_properties (bppTools PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif (BUILD_STATIC)

foreach (target ${bppsuite-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC} bppTools Eigen3::Eigen)
    set_target_properties (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED} bppTools Eigen3::Eigen)
    set_target_properties (${target} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${bppsuite-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
