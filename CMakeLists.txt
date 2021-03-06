#   Copyright 2015 Patrick Putnam
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
cmake_minimum_required (VERSION 2.8)
PROJECT( CLOTHO C CXX )

SET (${PROJECT_NAME}_LIB_CORENAME "clotho" )

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

SET (${PROJECT_NAME}_VER_MAJ 0)
SET (${PROJECT_NAME}_VER_MIN 1)
SET (${PROJECT_NAME}_VER_RELEASE 0 )

SET (${PROJECT_NAME}_RESOURCES_DIR ${${PROJECT_NAME}_SOURCE_DIR}/cmake)
SET (${PROJECT_NAME}_SRC_DIR       ${${PROJECT_NAME}_SOURCE_DIR}/src)
SET (${PROJECT_NAME}_EXAMPLE_DIR       ${${PROJECT_NAME}_SOURCE_DIR}/examples)
SET( ${PROJECT_NAME}_INCLUDE_DIR   ${${PROJECT_NAME}_SOURCE_DIR}/include )
SET( ${PROJECT_NAME}_UNITTEST_DIR   ${${PROJECT_NAME}_SOURCE_DIR}/unittest )

SET (CMAKE_INSTALL_PREFIX "./clotho" CACHE PATH "Install path prefix, prepended onto install directories")

IF (NOT ${PROJECT_NAME}_INSTALL_BIN_DIR )
    SET (${PROJECT_NAME}_INSTALL_BIN_DIR bin )
ENDIF (NOT ${PROJECT_NAME}_INSTALL_BIN_DIR )

IF (NOT ${PROJECT_NAME}_INSTALL_LIB_DIR )
    SET (${PROJECT_NAME}_INSTALL_LIB_DIR lib )
ENDIF (NOT ${PROJECT_NAME}_INSTALL_LIB_DIR )

IF (NOT ${PROJECT_NAME}_INSTALL_INCLUDE_DIR )
    SET(${PROJECT_NAME}_INSTALL_INCLUDE_DIR include)
ENDIF (NOT ${PROJECT_NAME}_INSTALL_INCLUDE_DIR )

################################################################################
SET ( ${PROJECT_NAME}_PACKAGE "clotho" )
SET ( ${PROJECT_NAME}_PACKAGE_NAME "${PROJECT_NAME}" )
SET (${PROJECT_NAME}_PACKAGE_VERSION "${${PROJECT_NAME}_VER_MAJ}.${${PROJECT_NAME}_VER_MIN}.${${PROJECT_NAME}_VER_RELEASE}")
SET (${PROJECT_NAME}_PACKAGE_VER_MAJOR "${${PROJECT_NAME}_VER_MAJ}.${${PROJECT_NAME}_VER_MIN}")
SET (${PROJECT_NAME}_PACKAGE_VER_MINOR "${${PROJECT_NAME}_VER_RELEASE}")
SET (${PROJECT_NAME}_PACKAGE_STRING "${${PROJECT_NAME}_PACKAGE_NAME} ${${PROJECT_NAME}_PACKAGE_VERSION}")
SET (${PROJECT_NAME}_PACKAGE_BUGREPORT "putnampp@gmail.com")

MESSAGE( STATUS ${${PROJECT_NAME}_PACKAGE_STRING} )

IF (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED )
    SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY 
        ${PROJECT_BINARY_DIR}/bin CACHE PATH "Executable directory")
    SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY
        ${PROJECT_BINARY_DIR}/lib CACHE PATH "Libraries directory")
    SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY
        ${PROJECT_BINARY_DIR}/lib CACHE PATH "Static Libraries directory")
ELSE (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED )
    IF (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY )
        SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH} )
    ENDIF (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY )
ENDIF (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED )

######################################################################################
## PROJECT COMPILE OPTIONS
######################################################################################

OPTION(WITH_CUDA_EXAMPLES "Build the CUDA example programs" OFF )
OPTION(WITH_GSL "Build with GSL" OFF )
OPTION(WITH_TC_MALLOC "Build with TC Malloc" ON )
OPTION(WITH_PROFILER "Build with Gperftools profiler" ON )
OPTION(WITH_RT "Build with POSIX Realtime Extensions library" OFF )
OPTION(WITH_PTHREAD "Build with POSIX Thread library" OFF )

################################################################################

LIST(APPEND boost_comps program_options)
LIST(APPEND boost_comps filesystem)
LIST(APPEND boost_comps system)
LIST(APPEND boost_comps chrono)

LIST(APPEND boost_comps unit_test_framework)

LIST(APPEND boost_comps python)
LIST(APPEND boost_comps thread)
LIST(APPEND boost_comps log)


FIND_PACKAGE( Boost COMPONENTS ${boost_comps} REQUIRED )
IF( Boost_FOUND )
    INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
    MESSAGE (STATUS "Boost Include: " ${Boost_INCLUDE_DIR} )
    MESSAGE (STATUS "Boost Library: " ${Boost_LIBRARIES} )
ENDIF( Boost_FOUND )
#LIST( APPEND COMMON_LIBRARIES ${Boost_LIBRARIES} )

#MESSAGE (STATUS "Program options component: " ${Boost_PROGRAM_OPTIONS_LIBRARY})
#
#foreach( ITEM ${Boost_LIBRARIES} )
#    if( NOT (ITEM MATCHES "(.*)python(.*)") )
#        LIST(APPEND COMMON_LIBRARIES ${ITEM})
#    endif()
#endforeach()

FIND_PACKAGE( PythonLibs 2.7 REQUIRED )
IF( PYTHONLIBS_FOUND )
    INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )
    MESSAGE( STATUS "Python Include: " ${PYTHON_INCLUDE_DIRS} )
ENDIF( PYTHONLIBS_FOUND )
#LIST( APPEND COMMON_LIBRARIES ${PYTHON_LIBRARIES})

IF( WITH_GSL )     
    FIND_PACKAGE( GSL )
    IF(GSL_FOUND )
        MESSAGE( STATUS "GSL Library: " ${GSL_LIBRARIES} )
        LIST( APPEND COMMON_LIBRARIES ${GSL_LIBRARIES} )
    ENDIF( GSL_FOUND )
ENDIF( WITH_GSL )

IF( WITH_CUDA_EXAMPLES )
    FIND_PACKAGE( CUDA REQUIRED )
    MESSAGE( STATUS "CXX: " ${CMAKE_CXX_COMPILER})
    SET(CUDA_HOST_COMPILER g++)
    SET(CUDA_SEPARABLE_COMPILATION ON )
    SET(CUDA_PROPAGATE_HOST_FLAGS OFF)
    #SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-gencode=arch=compute_30,code=sm_30")
    LIST(APPEND CUDA_NVCC_FLAGS "-g;-O3;-gencode=arch=compute_52,code=sm_52;-std=c++11;-lineinfo;-src-in-ptx")

    MESSAGE( STATUS "CUDA VERSION: " ${CUDA_VERSION} )
    MESSAGE( STATUS "CUDA LIBRARIES: " ${CUDA_LIBRARIES} )
    MESSAGE( STATUS "CUDA LIBRARIES: " ${CUDA_curand_LIBRARY} )
    MESSAGE( STATUS "CUDA 64 bit device code: " ${CUDA_64_BIT_DEVICE_CODE} )
    MESSAGE( STATUS "CUDA NVCC FLAGS: " ${CUDA_NVCC_FLAGS} )
ENDIF( WITH_CUDA_EXAMPLES )

SET( Gperftools_DIR "/opt/gperftools/latest" )

IF( WITH_TC_MALLOC )
    find_library( TCMALLOC_LIBRARY NAMES tcmalloc PATHS ${Gperftools_DIR}/lib )
    MESSAGE( STATUS "TCMALLOC: " ${TCMALLOC_LIBRARY})
    LIST( APPEND COMMON_LIBRARIES ${TCMALLOC_LIBRARY} )
ENDIF( WITH_TC_MALLOC )

IF( WITH_PROFILER )
    find_library( PROFILER_LIBRARY NAMES profiler PATHS ${Gperftools_DIR}/lib )
    MESSAGE( STATUS "PROFILER: " ${PROFILER_LIBRARY})
    LIST( APPEND COMMON_LIBRARIES ${PROFILER_LIBRARY} )
ENDIF( WITH_PROFILER )

IF( WITH_RT )
    find_library( REALTIME_LIB NAMES rt )
    LIST( APPEND COMMON_LIBRARIES ${REALTIME_LIB} )
ENDIF( WITH_RT )

IF( WITH_PTHREAD )
    FIND_PACKAGE( Threads )
    IF( NOT Threads_FOUND )
        MESSAGE(STATUS "Threads not found; Attempting to locate pthread library manually" )   
        find_library( PTHREAD_LIB NAMES pthread )
        IF( PTHREAD_LIB_NOTFOUND )
            MESSAGE( STATUS "Pthread not found" )
        ENDIF( PTHREAD_LIB_NOTFOUND )
    ENDIF( NOT Threads_FOUND )
ENDIF( WITH_PTHREAD )

FOREACH( L ${COMMON_LIBRARIES} ) 
    MESSAGE( STATUS ${L} )
ENDFOREACH(L)


################################################################################

#OPTION( LOG       "Compile with debugging; sets LOG_LEVEL to the default logging level"        OFF )
#OPTION( LOG_LEVEL "Set the default debug level; (-1 = OFF, 0 - Critical Message, 1 - Errors, 2 - Warnings, 3 - Status (default), 4 - Notes, 5 - Debug, 6 - All"   0 )  # OFF

#SET(MAX_ALLELES 256 CACHE INTEGER "Maximum alleles per locus")
#IF( MAX_ALLELES LESS 0 )
#    MESSAGE(FATAL_ERROR "Maximum alleles cannot be negative" )
#ENDIF()
#

SET(PWORD_SIZE 64 CACHE INTEGER "Bits per processor word")
IF(PWORD_SIZE LESS 0 )
    MESSAGE(FATAL_ERROR "Processor word size cannot be negative")
ENDIF()

IF( EXISTS "${${PROJECT_NAME}_INCLUDE_DIR}/config.h" )
    UNSET( ${${PROJECT_NAME}_INCLUDE_DIR}/config.h CACHE )
ENDIF()

#CONFIGURE_FILE( ${${PROJECT_NAME}_INCLUDE_DIR}/config.h.in ${${PROJECT_NAME}_INCLUDE_DIR}/config.h)
#CONFIGURE_FILE( ${${PROJECT_NAME}_INCLUDE_DIR}/cchmc/config.hpp.in ${${PROJECT_NAME}_INCLUDE_DIR}/cchmc/config.hpp)

#CONFIGURE_FILE( ${${PROJECT_NAME}_UNITTEST_DIR}/test_main.cpp.in ${${PROJECT_NAME}_UNITTEST_DIR}/test_main.cpp )

FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile @ONLY)
ADD_CUSTOM_TARGET(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/docs/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
ENDIF(DOXYGEN_FOUND)

set(CMAKE_CXX_FLAGS "-g -fPIC -Wall -O3 -std=c++11 -march=native -pedantic" )

LINK_DIRECTORIES( ${PROJECT_BINARY_DIR}/bin )

add_subdirectory( src )
add_subdirectory( examples )
add_subdirectory( unittest )
