# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 22/08/2009

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bppsuite CXX)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()

SET(CMAKE_CXX_FLAGS "-Wall -Weffc++ -Wshadow -Wconversion")
IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF()

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
ENDIF()

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NOT DOC_COMPRESS)
  SET(DOC_COMPRESS gzip CACHE STRING
      "Set program for compressing documentation."
      FORCE)
ENDIF(NOT DOC_COMPRESS)

IF(NOT DOC_COMPRESS_EXT)
  SET(DOC_COMPRESS_EXT gz CACHE STRING
      "Set extension of compressed documentation."
      FORCE)
ENDIF(NOT DOC_COMPRESS_EXT)


IF(NO_DEP_CHECK)
  MESSAGE("-- Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)


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
IF(MAN)
  FIND_PROGRAM(NROFF_EXE NAMES nroff)
  IF(NROFF_EXE)
    MESSAGE("-- Found nroff here: ${NROFF_EXE}")
    MESSAGE("   Adding targets: man")

    ADD_CUSTOM_TARGET(man
      ALL
      COMMAND cp bppml.1.txt bppml.1
      COMMAND ${DOC_COMPRESS} -f bppml.1
      COMMAND cp bppseqgen.1.txt bppseqgen.1
      COMMAND ${DOC_COMPRESS} -f bppseqgen.1
      COMMAND cp bppdist.1.txt bppdist.1
      COMMAND ${DOC_COMPRESS} -f bppdist.1
      COMMAND cp bpppars.1.txt bpppars.1
      COMMAND ${DOC_COMPRESS} -f bpppars.1
      COMMAND cp bppseqman.1.txt bppseqman.1
      COMMAND ${DOC_COMPRESS} -f bppseqman.1
      COMMAND cp bppconsense.1.txt bppconsense.1
      COMMAND ${DOC_COMPRESS} -f bppconsense.1
      COMMAND cp bppancestor.1.txt bppancestor.1
      COMMAND ${DOC_COMPRESS} -f bppancestor.1
      COMMAND cp bppreroot.1.txt bppreroot.1
      COMMAND ${DOC_COMPRESS} -f bppreroot.1
      COMMAND cp bppphysamp.1.txt bppphysamp.1
      COMMAND ${DOC_COMPRESS} -f bppphysamp.1
      COMMAND cp bpptreedraw.1.txt bpptreedraw.1
      COMMAND ${DOC_COMPRESS} -f bpptreedraw.1
      COMMAND cp bppalnscore.1.txt bppalnscore.1
      COMMAND ${DOC_COMPRESS} -f bppalnscore.1
      COMMAND cp bppmixedlikelihoods.1.txt bppmixedlikelihoods.1
      COMMAND ${DOC_COMPRESS} -f bppmixedlikelihoods.1
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
      )
    SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "man/bppml.1.${DOC_COMPRESS_EXT};man/bppseqgen.1.${DOC_COMPRESS_EXT};man/bppdist.1.${DOC_COMPRESS_EXT};man/bpppars.1.${DOC_COMPRESS_EXT};man/bppseqman.1.${DOC_COMPRESS_EXT};man/bppconsense.1.${DOC_COMPRESS_EXT};man/bppancestor.1.${DOC_COMPRESS_EXT};man/bppreroot.1.${DOC_COMPRESS_EXT};man/bppphysamp.1.${DOC_COMPRESS_EXT};man/bpptreedraw.1.${DOC_COMPRESS_EXT};man/bppalnscore.1.${DOC_COMPRESS_EXT};man/bppmixedlikelihoods.1.${DOC_COMPRESS_EXT}")

  ELSE()
    MESSAGE(FATAL_ERROR "Program nroff required but not found.")
  ENDIF()
ENDIF(MAN)

IF(INFO)
  FIND_PROGRAM(MAKEINFO_EXE NAMES makeinfo)
  IF(MAKEINFO_EXE)
    MESSAGE("-- Found makeinfo here: ${MAKEINFO_EXE}")
    MESSAGE("   Adding targets: info, html")
 
    SET(ADD_INFO_TO "ALL")
    MESSAGE("   Adding target info to target all")
    
    ADD_CUSTOM_TARGET(info
      ${ADD_INFO_TO}
      COMMAND ${MAKEINFO_EXE} bppsuite.texi
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
      )

    ADD_CUSTOM_TARGET(html
      COMMAND ${MAKEINFO_EXE} --html --css-ref=http://www.w3.org/StyleSheets/Core/Steely bppsuite.texi
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
      )
  ELSE(MAKEINFO_EXE)
    MESSAGE(FATAL_ERROR p"Program makeinfo required but not found.")
  ENDIF(MAKEINFO_EXE)
ENDIF(INFO)

FIND_PROGRAM(PDFTEX_EXE NAMES pdftex)
IF(PDFTEX_EXE)
  MESSAGE("-- Found pdftex here: ${PDFTEX_EXE}")
  MESSAGE("   Adding target: pdf")

  ADD_CUSTOM_TARGET(pdf
    COMMAND ${PDFTEX_EXE} bppsuite.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "doc/bppsuite.pdf;doc/bppsuite.aux;doc/bppsuite.cp;doc/bppsuite.fn;doc/bppsuite.info;doc/bppsuite.ky;doc/bppsuite.log;doc/bppsuite.pg;doc/bppsuite.toc;doc/bppsuite.tp;doc/bppsuite.vr")
ENDIF(PDFTEX_EXE)

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
IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Bpp/Phyl/Tree/Tree.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-core Bpp/Clonable.h)

# Subdirectories
ADD_SUBDIRECTORY(bppSuite)
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(man)

ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "bppsuite")
SET(CPACK_PACKAGE_VENDOR "Bio++ Development Team")
SET(CPACK_PACKAGE_VERSION "2.2.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "2")
SET(CPACK_PACKAGE_VERSION_MINOR "2")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
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
 ".*\\\\.dmg"
 ".*\\\\.sh"
 ".*\\\\..*\\\\.swp"
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
 "bppSuite/bppalnscore"
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
 "man/.*\\\\.1.${DOC_COMPRESS_EXT}"
 "debian/tmp"
 "debian/bppsuite/"
 "debian/bppsuite\\\\.substvars"
 "debian/bppsuite\\\\.debhelper"
 "debian/debhelper\\\\.log"
 "install_manifest.txt"
 "DartConfiguration.tcl"
 ${CPACK_SOURCE_IGNORE_FILES}
)
IF (MACOS)
  SET(CPACK_GENERATOR "Bundle")
ENDIF()

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
# 'clean' is not (yet) a first class target. However, we need to clean the directories before building the sources:
IF("${CMAKE_GENERATOR}" MATCHES "Make")
  ADD_CUSTOM_TARGET(make_clean
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
  )
  ADD_DEPENDENCIES(dist make_clean)

  ADD_CUSTOM_TARGET(make_clean_man
    COMMAND rm -f *.${DOC_COMPRESS_EXT}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/man
  )
  ADD_DEPENDENCIES(dist make_clean_man)
ENDIF()

IF(NOT NO_DEP_CHECK)
IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(origdist COMMAND cp ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz ../${CPACK_DEBSOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(origdist dist)
ADD_CUSTOM_TARGET(deb dpkg-buildpackage -uc -us -i${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(deb origdist)
ADD_DEPENDENCIES(deb info man)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info man)

#ADD_DEPENDENCIES(info install)

ENDIF(UNIX)

ENDIF(NOT NO_DEP_CHECK)
