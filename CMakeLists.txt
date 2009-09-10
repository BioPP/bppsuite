# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 22/08/2009

# Global parameters
cmake_minimum_required(VERSION 2.6)
project(bppsuite CXX)
set(CMAKE_BUILD_TYPE "Release")

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
if(CMAKE_INSTALL_PREFIX)
  set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif(CMAKE_INSTALL_PREFIX)

# Find the bpp-utils library
find_path(BPPUTILS_INCLUDE_DIR Utils/Clonable.h)
set(BPPUTILS_NAMES bpp-utils bpp-utilslib bpp-utilsdll)
find_library(BPPUTILS_LIBRARY NAMES ${BPPUTILS_NAMES})
set( BPPUTILS_LIBRARIES ${BPPUTILS_LIBRARY} )
if(BPPUTILS_LIBRARY)
  message("bpp-utils found here:")
  message("  includes : ${BPPUTILS_INCLUDE_DIR}")
  message("  libraries: ${BPPUTILS_LIBRARIES}")
else(BPPUTILS_LIBRARY)
  message(FATAL_ERROR "bpp-utils required but not found")
endif(BPPUTILS_LIBRARY)

# Include the bpp-utils dependency
include_directories(${BPPUTILS_INCLUDE_DIR})
set(LIBS ${LIBS} ${BPPUTILS_LIBRARIES})

# Find the bpp-numcalc library
find_path(BPPNUMCALC_INCLUDE_DIR NumCalc/NumTools.h)
set(BPPNUMCALC_NAMES bpp-numcalc bpp-numcalclib bpp-numcalcdll)
find_library(BPPNUMCALC_LIBRARY NAMES ${BPPNUMCALC_NAMES})
set( BPPNUMCALC_LIBRARIES ${BPPNUMCALC_LIBRARY} )
if(BPPNUMCALC_LIBRARY)
  message("bpp-numcalc found here:")
  message("  includes : ${BPPNUMCALC_INCLUDE_DIR}")
  message("  libraries: ${BPPNUMCALC_LIBRARIES}")
else(BPPNUMCALC_LIBRARY)
  message(FATAL_ERROR "bpp-numcalc required but not found")
endif(BPPNUMCALC_LIBRARY)

# Include the bpp-numcalc dependency
include_directories(${BPPNUMCALC_INCLUDE_DIR})
set(LIBS ${LIBS} ${BPPNUMCALC_LIBRARIES})

# Find the bpp-seq library
find_path(BPPSEQ_INCLUDE_DIR Seq/Alphabet.h)
set(BPPSEQ_NAMES bpp-seq bpp-seqlib bpp-seqdll)
find_library(BPPSEQ_LIBRARY NAMES ${BPPSEQ_NAMES})
set( BPPSEQ_LIBRARIES ${BPPSEQ_LIBRARY} )
if(BPPSEQ_LIBRARY)
  message("bpp-seq found here:")
  message("  includes : ${BPPSEQ_INCLUDE_DIR}")
  message("  libraries: ${BPPSEQ_LIBRARIES}")
else(BPPSEQ_LIBRARY)
  message(FATAL_ERROR "bpp-seq required but not found")
endif(BPPSEQ_LIBRARY)

# Include the bpp-seq dependency
include_directories(${BPPSEQ_INCLUDE_DIR})
set(LIBS ${LIBS} ${BPPSEQ_LIBRARIES})

# Find the bpp-phyl library
find_path(BPPPHYL_INCLUDE_DIR Phyl/Tree.h)
set(BPPPHYL_NAMES bpp-phyl bpp-phyllib bpp-phyldll)
find_library(BPPPHYL_LIBRARY NAMES ${BPPPHYL_NAMES})
set( BPPPHYL_LIBRARIES ${BPPPHYL_LIBRARY} )
if(BPPPHYL_LIBRARY)
  message("bpp-phyl found here:")
  message("  includes : ${BPPPHYL_INCLUDE_DIR}")
  message("  libraries: ${BPPPHYL_LIBRARIES}")
else(BPPPHYL_LIBRARY)
  message(FATAL_ERROR "bpp-phyl required but not found")
endif(BPPPHYL_LIBRARY)

# Include the bpp-phyl dependency
include_directories(${BPPPHYL_INCLUDE_DIR})
set(LIBS ${LIBS} ${BPPPHYL_LIBRARIES})

# We should check for the Texinfo dependecy, but no FindPackage is available here...
add_custom_target(info
  COMMAND makeinfo bppsuite.texi
  WORKING_DIRECTORY doc
  )
add_custom_target(html
  COMMAND makeinfo --html --css-ref=http://www.w3.org/StyleSheets/Core/Steely bppsuite.texi
  WORKING_DIRECTORY doc
  )

# Subdirectories
add_subdirectory(bppSuite)
add_subdirectory(doc)

# Packager
include(CPack)
