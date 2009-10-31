# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 22/08/2009

# Global parameters
cmake_minimum_required(VERSION 2.6)
project(bppsuite CXX)
if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT DEFINED CMAKE_BUILD_TYPE)  
set(CMAKE_CXX_FLAGS "-Wall")

#static linkage?
if(DEFINED BUILD_STATIC)
  if(BUILD_STATIC)
    message("-- Static linkage requested.")
    set(CMAKE_CXX_FLAGS "-static ${CMAKE_CXX_FLAGS}")
  endif(BUILD_STATIC)
else(DEFINED BUILD_STATIC)
  set(BUILD_STATIC 0) #Dynamic linkage by default
endif(DEFINED BUILD_STATIC)

#here is a useful function:
macro(improved_find_library OUTPUT_LIBS lib_name include_to_find static_lib)
  #start:
  find_path(${lib_name}_INCLUDE_DIR ${include_to_find})

  if(${static_lib})
    find_path(${lib_name}_LIBRARY_DIR lib${lib_name}.a HINTS ${CMAKE_PREFIX_PATH} ${CMAKE_PREFIX_PATH}/lib)
    if(${lib_name}_LIBRARY_DIR)
      set(${lib_name}_LIBRARY ${${lib_name}_LIBRARY_DIR}/lib${lib_name}.a)
    else(${lib_name}_LIBRARY_DIR)
      message(FATAL_ERROR "${lib_name} static library could not be found.")
    endif(${lib_name}_LIBRARY_DIR)
  else(${static_lib})
    set(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
    find_library(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})
  endif(${static_lib})

  if(${lib_name}_LIBRARY)
    message("-- Library ${lib_name} found here:")
    message("   includes : ${${lib_name}_INCLUDE_DIR}")
    message("   libraries: ${${lib_name}_LIBRARY}")
  else(${lib_name}_LIBRARY)
    message(FATAL_ERROR "${lib_name} required but not found.")
  endif(${lib_name}_LIBRARY)
  
  #add the dependency:
  include_directories(${${lib_name}_INCLUDE_DIR})
  set(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
endmacro(improved_find_library)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
if(CMAKE_INSTALL_PREFIX)
  set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif(CMAKE_INSTALL_PREFIX)

#Find the libraries. The order is very important for staitc linkage, it won't
#link if you change it!
improved_find_library(LIBS bpp-phyl Phyl/Tree.h ${BUILD_STATIC})
improved_find_library(LIBS bpp-seq Seq/Alphabet.h ${BUILD_STATIC})
improved_find_library(LIBS bpp-numcalc NumCalc/NumTools.h ${BUILD_STATIC})
improved_find_library(LIBS bpp-utils Utils/Clonable.h ${BUILD_STATIC})

# We should check for the Texinfo dependecy, but no FindPackage is available here...
add_custom_target(info
  COMMAND makeinfo bppsuite.texi
  WORKING_DIRECTORY doc
  )
add_custom_target(pdf
  COMMAND pdftex bppsuite.texi
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
