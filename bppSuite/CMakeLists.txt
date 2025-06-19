# SPDX-FileCopyrightText: The Bio++ Development Group
#
# SPDX-License-Identifier: CECILL-2.1


# Executables of bppsuite.
# Generation of targets from file name is not automated in case of executables not following the pattern.

add_executable (bppml bppML.cpp)
add_executable (bppseqgen bppSeqGen.cpp)
add_executable (bppdist bppDist.cpp)
add_executable (bpppars bppPars.cpp)
add_executable (bppseqman bppSeqMan.cpp)
add_executable (bppconsense bppConsense.cpp)
add_executable (bppancestor bppAncestor.cpp)
add_executable (bppmixedlikelihoods bppMixedLikelihoods.cpp)
add_executable (bppbranchlik bppBranchLik.cpp)
add_executable (bppreroot bppReRoot.cpp)
add_executable (bpptreedraw bppTreeDraw.cpp)
add_executable (bppalnscore bppAlnScore.cpp)
add_executable (bpppopstats bppPopStats.cpp)
add_executable (bppr4s bppR4S.cpp)

set (bppsuite-targets
  bppml
  bppseqgen
  bppdist
  bpppars
  bppseqman
  bppconsense
  bppancestor
  bppmixedlikelihoods
  bppbranchlik
  bppreroot
  bpptreedraw
  bppalnscore
  bpppopstats
  bppr4s
  )

foreach (target ${bppsuite-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC} Eigen3::Eigen)
    set_target_properties (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED} Eigen3::Eigen)
    set_target_properties (${target} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${bppsuite-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
