# SPDX-FileCopyrightText: The Bio++ Development Group
#
# SPDX-License-Identifier: CECILL-2.1

# CMake script for Bio++ Program Suite
# Authors:
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 22/08/2009

cmake_minimum_required (VERSION 3.10)
project (bppsuite CXX)

# Compile options
add_compile_options(-std=c++14 -Wall -Weffc++ -Wshadow -Wconversion)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NO_DEP_CHECK)
  MESSAGE(STATUS "Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz DEPENDS dist info)

#static linkage?
if(NOT BUILD_STATIC)
  set(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
endif()
if(BUILD_STATIC)
  message(STATUS "Static linkage requested.")
  add_compile_options(-static -static-libgcc)
endif()

# Check compression program
# COMPRESS_PROGRAM controls the choice of program
# COMPRESS_EXT can be used to override the file extension
if (NOT COMPRESS_PROGRAM)
  set (COMPRESS_PROGRAM gzip CACHE STRING "Set program for compressing documentation" FORCE)
endif ()
find_program (COMPRESS_BIN NAMES ${COMPRESS_PROGRAM} DOC "${COMPRESS_PROGRAM} compression program")
if (NOT COMPRESS_BIN)
  message (STATUS "${COMPRESS_PROGRAM} program not found, text doc will not be compressed")
else ()
  # Deduce COMPRESS_EXT for known compression programs if not set
  if (NOT COMPRESS_EXT)
    if (${COMPRESS_PROGRAM} STREQUAL "gzip")
      set (COMPRESS_EXT "gz")
    elseif (${COMPRESS_PROGRAM} STREQUAL "bzip2")
      set (COMPRESS_EXT "bz2")
    else ()
      set (COMPRESS_EXT "${COMPRESS_PROGRAM}") # Default: program name (works for xz/lzma)
    endif ()
  endif ()
  # Generate command line args (always add -c to output compressed file to stdout)
  if (${COMPRESS_PROGRAM} STREQUAL "gzip")
    # -n for no timestamp in files (reproducible builds)
    # -9 for maximum compression (lintian error)
    set (COMPRESS_ARGS -c -n -9)
  else ()
    set (COMPRESS_ARGS -c)
  endif ()
  message (STATUS "Found ${COMPRESS_BIN} compression program, using file extension .${COMPRESS_EXT}")
endif ()

# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

include (GNUInstallDirs)
find_package (bpp-phyl3 1.0.0 REQUIRED)
find_package (bpp-popgen3 1.0.0 REQUIRED)

find_package (Eigen3 3.3 REQUIRED PATHS /usr/lib NO_MODULE)

# Subdirectories
add_subdirectory (bppSuite)
add_subdirectory (doc)

ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "bppsuite")
SET(CPACK_PACKAGE_VENDOR "Bio++ Development Team")
SET(CPACK_PACKAGE_VERSION "3.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "3")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Bio++ Program Suite")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSES/CECILL-2.1.txt")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
SET(CPACK_SOURCE_IGNORE_FILES 
       "/build/" 
       "/\\\\.git/" 
       "/\\\\.gitignore" 
       ${CPACK_SOURCE_IGNORE_FILES}
       )

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

