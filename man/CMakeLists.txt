# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 22/08/2009

IF(MAN)
  INSTALL(FILES bppml.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppseqgen.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppancestor.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bpppars.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppdist.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppconsense.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppseqman.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppreroot.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppphysamp.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bpptreedraw.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppalnscore.1.gz DESTINATION share/man/man1)
  INSTALL(FILES bppmixedlikelihoods.1.gz DESTINATION share/man/man1)
ENDIF(MAN)
