

# DBC Note these are hardwared to the local install directory
set(PALISADE_INCLUDE "${CMAKE_INSTALL_PREFIX}/include/palisade")
set(PALISADE_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")
message(DEBUG "PALISADE_INCLUDE: " ${PALISADE_INCLUDE})
message(DEBUG "PALISADE_LIB: " ${PALISADE_INCLUDE})

include_directories( ${PALISADE_INCLUDE}/core)
include_directories( ${PALISADE_INCLUDE}/pke)
include_directories( ${PALISADE_INCLUDE}/subgaussian)
include_directories( include )
include_directories( lib )

include_directories( ../src/trapdoor/include )
include_directories( ../src/trapdoor/lib )
include_directories( ../src/abe/include )
include_directories( ../src/abe/lib )
include_directories( ../src/signature/include )
include_directories( ../src/signature/lib )

set( REGEX_FLAG -DHAVE_STD_REGEX )

if ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" )
	set( WARNING_FLAG "" )
else()
	set( WARNING_FLAG -Wno-unused-but-set-variable )
endif()

set (BMLIBS PUBLIC PALISADEtrapdoor PUBLIC PALISADEabe PUBLIC PALISADEsignature PUBLIC PALISADEsubgaussian PALISADEpke PUBLIC PALISADEbinfhe PUBLIC PALISADEcore ${THIRDPARTYLIBS} PUBLIC benchmark ${OpenMP_CXX_FLAGS})

set (BMAPPS "")
file (GLOB BMARK_SRC_FILES CONFIGURE_DEPENDS src/*.cpp)
foreach (app ${BMARK_SRC_FILES})
	get_filename_component ( exe ${app} NAME_WE )
	add_executable ( ${exe} ${app} )
	set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmark)
	target_include_directories( ${exe} PUBLIC ../third-party/google-benchmark/include )
	set_target_properties( ${exe} PROPERTIES COMPILE_FLAGS
		"${REGEX_FLAG} -DHAVE_STEADY_CLOCK -DNDEBUG ${WARNING_FLAG} -Wno-unused-variable" )
	target_link_libraries ( ${exe} ${BMLIBS} )
	set (BMAPPS ${BMAPPS} ${exe})
endforeach()

add_custom_target( allbenchmark )
if (${BMAPPS})
  add_dependencies( allbenchmark ${BMAPPS} )
endif()

add_custom_command( OUTPUT benchmarkinfocmd COMMAND echo Builds benchmark lib and these apps: ${BMAPPS} )
add_custom_target( benchmarkinfo DEPENDS benchmarkinfocmd )
