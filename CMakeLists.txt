# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 22/08/2009

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)
PROJECT(bppsuite CXX)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()

SET(CMAKE_CXX_FLAGS "-Wall -Weffc++ -Wshadow")
IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF()

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
ENDIF()

#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE("-- Static linkage requested.")
  SET(CMAKE_CXX_FLAGS "-static -static-libgcc ${CMAKE_CXX_FLAGS}")
ENDIF()

#build info?
IF(NOT DEFINED INFO)
  SET(INFO UNIX)
ENDIF(NOT DEFINED INFO)

#build man pages?
IF(NOT DEFINED MAN)
  SET(MAN UNIX)
ENDIF(NOT DEFINED MAN)

#find executables for documentation
FIND_PROGRAM(MAKEINFO_EXE NAMES makeinfo)
FIND_PROGRAM(PDFTEX_EXE NAMES pdftex)
FIND_PROGRAM(NROFF_EXE NAMES nroff)

#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})

  IF(BUILD_STATIC)
    SET(${lib_name}_STATIC_NAMES lib${lib_name}.a)
    FIND_LIBRARY(${lib_name}_STATIC_LIBRARY NAMES ${${lib_name}_STATIC_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_STATIC_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   static libraries: ${${lib_name}_STATIC_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_STATIC_LIBRARY})
  ELSE()
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_LIBRARY)
      MESSAGE("-- Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
  ENDIF()
  
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

#Find the libraries. The order is very important for static linkage, it won't
#link if you change it!
IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Bpp/Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-core Bpp/Clonable.h)

IF(MAKEINFO_EXE)
  MESSAGE("-- Found makeinfo here: ${MAKEINFO_EXE}")
  MESSAGE("   Adding targets: info, html")

  IF(INFO)
    SET(ADD_INFO_TO "ALL")
    MESSAGE("   Adding target info to target all")
  ENDIF(INFO)

  ADD_CUSTOM_TARGET(info
    ${ADD_INFO_TO}
    COMMAND ${MAKEINFO_EXE} bppsuite.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )

  ADD_CUSTOM_TARGET(html
    COMMAND ${MAKEINFO_EXE} --html --css-ref=http://www.w3.org/StyleSheets/Core/Steely bppsuite.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
ENDIF(MAKEINFO_EXE)

IF(PDFTEX_EXE)
  MESSAGE("-- Found pdftex here: ${PDFTEX_EXE}")
  MESSAGE("   Adding target: pdf")

  ADD_CUSTOM_TARGET(pdf
    COMMAND ${PDFTEX_EXE} bppsuite.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
ENDIF(PDFTEX_EXE)

IF(NROFF_EXE)
  MESSAGE("-- Found nroff here: ${NROFF_EXE}")
  MESSAGE("   Adding targets: man, man-gz")

  ADD_CUSTOM_TARGET(man
    COMMAND cp bppml.1.txt bppml.1
    COMMAND cp bppseqgen.1.txt bppseqgen.1
    COMMAND cp bppdist.1.txt bppdist.1
    COMMAND cp bpppars.1.txt bpppars.1
    COMMAND cp bppseqman.1.txt bppseqman.1
    COMMAND cp bppconsense.1.txt bppconsense.1
    COMMAND cp bppancestor.1.txt bppancestor.1
    COMMAND cp bppreroot.1.txt bppreroot.1
    COMMAND cp bppphysamp.1.txt bppphysamp.1
    COMMAND cp bpptreedraw.1.txt bpptreedraw.1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
    )

  ADD_CUSTOM_TARGET(man-gz ALL
    COMMAND gzip -f bppml.1
    COMMAND gzip -f bppseqgen.1
    COMMAND gzip -f bppdist.1
    COMMAND gzip -f bpppars.1
    COMMAND gzip -f bppseqman.1
    COMMAND gzip -f bppconsense.1
    COMMAND gzip -f bppancestor.1
    COMMAND gzip -f bppreroot.1
    COMMAND gzip -f bppphysamp.1
    COMMAND gzip -f bpptreedraw.1
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
    )
  
  ADD_DEPENDENCIES(man-gz man)
ENDIF(NROFF_EXE)

# Subdirectories
ADD_SUBDIRECTORY(bppSuite)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(man)

# Packager
SET(CPACK_PACKAGE_NAME "bppsuite")
SET(CPACK_PACKAGE_VENDOR "Bio++ Development Team")
SET(CPACK_PACKAGE_VERSION "0.6.1")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "6")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Bio++ Program Suite")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 ".*\\\\.svn"
 ".*\\\\.gz"
 ".*\\\\.zip"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 "\\\\..*\\\\.swp"
 ".*stamp"
 "\\\\.sh"
 "bppSuite/\\\\..*"
 "bppSuite/bppancestor"
 "bppSuite/bppml"
 "bppSuite/bppphysamp"
 "bppSuite/bppseqman"
 "bppSuite/bppconsense"
 "bppSuite/bppreroot"
 "bppSuite/bpppars"
 "bppSuite/bpptreedraw"
 "bppSuite/bppdist"
 "bppSuite/bppseqgen"
 "doc/bppsuite/"
 "doc/bppsuite\\\\.info"
 "doc/bppsuite\\\\.toc"
 "doc/bppsuite\\\\.vr"
 "doc/bppsuite\\\\.tp"
 "doc/bppsuite\\\\.log"
 "doc/bppsuite\\\\.fn"
 "doc/bppsuite\\\\.ky"
 "doc/bppsuite\\\\.cp"
 "doc/bppsuite\\\\.pg"
 "doc/bppsuite\\\\.aux"
 "doc/bppsuite\\\\.pdf"
 "man/.*\\\\.1.gz"
 "debian/bppsuite/"
 "debian/bppsuite\\\\.substvars"
 "debian/bppsuite\\\\.debhelper"
 "debian/debhelper\\\\.log"
 "debian/tmp"
 ${CPACK_SOURCE_IGNORE_FILES}
)
IF (MACOS)
  SET(CPACK_GENERATOR "Bundle")
ENDIF()

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(deb dpkg-buildpackage)
ADD_DEPENDENCIES(deb dist info man)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info man-gz)

#ADD_DEPENDENCIES(info install)

ENDIF(UNIX)

