cmake_minimum_required (VERSION 3.5.1)

# this version will get incremented only if new changes are made in this repo
set(PALISADE_VERSION_MAJOR 1)
set(PALISADE_VERSION_MINOR 11)
set(PALISADE_VERSION_PATCH 0)
set(PALISADE_VERSION ${PALISADE_VERSION_MAJOR}.${PALISADE_VERSION_MINOR}.${PALISADE_VERSION_PATCH})

project(EXAMPLES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#change to 1.11 once v1.11 is released
find_package(Palisade 1.10.6...)

set( CMAKE_CXX_FLAGS ${PALISADE_CXX_FLAGS} )
set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/third-party/include )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )

# third party directories

include_directories( third-party/google-test/googletest third-party/google-test/googletest/include )

### add directories for other PALISADE modules

link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
link_libraries( ${PALISADE_LIBRARIES} )
message(STATUS "PALISADE_LIBRARIES found: ${PALISADE_LIBRARIES}")

option( BUILD_UNITTESTS "Set to ON to build unit tests for the library" ON)
option( BUILD_EXAMPLES "Set to ON to build examples for the library" ON)

message( STATUS "BUILD_UNITTESTS:  ${BUILD_UNITTESTS}")
message( STATUS "BUILD_EXAMPLES:   ${BUILD_EXAMPLES}")

message("***** INSTALL IS AT ${CMAKE_INSTALL_PREFIX}; to change, run cmake with -DCMAKE_INSTALL_PREFIX=/your/path")
set (CMAKE_INSTALL_MESSAGE LAZY)

# when running "make clean", additionally deletes the CMake cache file
set(ADDITIONAL_CLEAN_FILES "")
LIST(APPEND ADDITIONAL_CLEAN_FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt)

# Add the additional "make clean" files
GET_DIRECTORY_PROPERTY(clean_files ADDITIONAL_CLEAN_FILES)
LIST(APPEND            clean_files ${ADDITIONAL_CLEAN_FILES})
LIST(REMOVE_DUPLICATES clean_files)
LIST(REMOVE_ITEM       clean_files "")
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_CLEAN_FILES "${clean_files}")
message(STATUS "clean files ${clean_files}")

add_subDirectory(src/abe)
