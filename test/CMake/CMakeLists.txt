cmake_minimum_required (VERSION 3.5.1)

# Create the abe project
project(PALISADEabeTest CXX)

# Use C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Palisade 1.11.0 REQUIRED)

set( CMAKE_CXX_FLAGS ${PALISADE_CXX_FLAGS} )
set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/third-party/include )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )

link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )

link_libraries( ${PALISADE_SHARED_LIBRARIES} )

find_package(PALISADEabe 1.11.0 REQUIRED)

if (PALISADEabe_FOUND)
  message(NOTICE "Libraries for mypackage found")
  message(NOTICE "Package Version: " ${PALISADEabe_VERSION})
else()
  message(FATAL_ERROR "Failed to find mypackage.")
endif()

include_directories( ${ABE_INCLUDE} )
link_directories( ${ABE_LIBDIR} )
link_libraries( ${ABE_SHARED_LIBRARIES} )

message(NOTICE "ABE_INCLUDE: " ${ABE_INCLUDE})
message(NOTICE "ABE_LIBRARIES: " ${ABE_LIBRARIES})
message(NOTICE "ABE_SHARED_LIBRARIES: " ${ABE_SHARED_LIBRARIES})

##########################################################
### now build all files in examples and put executables in bin/abe/examples
##########################################################
#for all src files in examples directory
file (GLOB ABE_EXAMPLES_SRC_FILES CONFIGURE_DEPENDS ../../src/abe/examples/*.cpp)
foreach (app ${ABE_EXAMPLES_SRC_FILES})
  get_filename_component ( exe ${app} NAME_WE )
  add_executable ( ${exe} ${app} )
  set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
endforeach()
