#
# CMakeLists.txt for TRAPDOOR library
#


# all files named *-impl.cpp are compiled to form the library
file (GLOB TRAPDOOR_SRC_FILES CONFIGURE_DEPENDS lib/*-impl.cpp lib/*/*-impl.cpp)

# DBC Note these are hardwared to the local install directory
set(PALISADE_INCLUDE "${CMAKE_INSTALL_PREFIX}/include/palisade")
set(PALISADE_LIBDIR "${CMAKE_INSTALL_PREFIX}/lib")
message(DEBUG "PALISADE_INCLUDE: " ${PALISADE_INCLUDE})
message(DEBUG "PALISADE_LIB: " ${PALISADE_INCLUDE})

include_directories( ${PALISADE_INCLUDE}/core)
include_directories( ${PALISADE_INCLUDE}/pke)
include_directories( include )
include_directories( lib )

set(TRAPDOOR_VERSION_MAJOR ${PALISADE_VERSION_MAJOR})
set(TRAPDOOR_VERSION_MINOR ${PALISADE_VERSION_MINOR})
set(TRAPDOOR_VERSION_PATCH ${PALISADE_VERSION_PATCH})
set(TRAPDOOR_VERSION ${TRAPDOOR_VERSION_MAJOR}.${TRAPDOOR_VERSION_MINOR}.${TRAPDOOR_VERSION_PATCH})

add_library(trapdoorobj OBJECT ${TRAPDOOR_SRC_FILES})
#DBC removed thenext line
#add_dependencies(trapdoorobj PALISADEpke)
set_property(TARGET trapdoorobj PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library (PALISADEtrapdoor SHARED $<TARGET_OBJECTS:trapdoorobj>)
set_property(TARGET PALISADEtrapdoor PROPERTY VERSION ${TRAPDOOR_VERSION})
set_property(TARGET PALISADEtrapdoor PROPERTY SOVERSION ${TRAPDOOR_VERSION_MAJOR})
set_property(TARGET PALISADEtrapdoor PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
install(TARGETS PALISADEtrapdoor
	EXPORT PalisadeTargets
	DESTINATION lib)
if( ${BUILD_STATIC} )
	add_library (PALISADEtrapdoor_static STATIC $<TARGET_OBJECTS:trapdoorobj>)
	set_property(TARGET PALISADEtrapdoor_static PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
	install(TARGETS PALISADEtrapdoor_static
		EXPORT PalisadeTargets
		DESTINATION lib)
endif()

install(DIRECTORY include/
	DESTINATION include/palisade/trapdoor)

set (TRAPDOORLIBS PUBLIC PALISADEtrapdoor PUBLIC PALISADEpke PUBLIC PALISADEcore ${THIRDPARTYLIBS} ${OpenMP_CXX_FLAGS})

target_link_libraries (PALISADEtrapdoor PUBLIC PALISADEpke PUBLIC PALISADEcore ${THIRDPARTYLIBS} ${OpenMP_CXX_FLAGS})
if( ${BUILD_STATIC} )
	target_link_libraries (PALISADEtrapdoor_static PUBLIC PALISADEpke_static PUBLIC PALISADEcore_static ${THIRDPARTYSTATICLIBS} ${OpenMP_CXX_FLAGS})
endif()

message(DEBUG "TRAPDOORLIBS: " ${TRAPDOORLIBS})
message(DEBUG "PALISADEcore: " ${PALISADEcore})

add_custom_target( alltrapdoor )
add_dependencies( alltrapdoor PALISADEtrapdoor )

if( BUILD_UNITTESTS )
	file (GLOB TRAPDOOR_TEST_SRC_FILES CONFIGURE_DEPENDS unittest/*.cpp)
	add_executable (trapdoor_tests ${TRAPDOOR_TEST_SRC_FILES} ${UNITTESTMAIN})
	set_property(TARGET trapdoor_tests PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/unittest)
	target_link_libraries ( trapdoor_tests ${TRAPDOORLIBS} )
	if (NOT ${USE_OpenMP})
		target_link_libraries ( trapdoor_tests PRIVATE Threads::Threads)
	endif()

	add_dependencies( alltrapdoor trapdoor_tests )

	add_custom_command( OUTPUT runtrapdoortests WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND ${CMAKE_BINARY_DIR}/unittest/trapdoor_tests )
	add_custom_target( testtrapdoor DEPENDS trapdoor_tests runtrapdoortests )
endif()

set (TRAPDOORAPPS "")
if (BUILD_EXAMPLES)
	file (GLOB TRAPDOOR_EXAMPLES_SRC_FILES CONFIGURE_DEPENDS examples/*.cpp)
	foreach (app ${TRAPDOOR_EXAMPLES_SRC_FILES})
		get_filename_component ( exe ${app} NAME_WE )
		add_executable ( ${exe} ${app} )
		set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/examples/trapdoor)
		set( TRAPDOORAPPS ${TRAPDOORAPPS} ${exe} )
		target_link_libraries ( ${exe} ${TRAPDOORLIBS} )
	endforeach()

	add_custom_target( alltrapdoorexamples )
	add_dependencies( alltrapdoorexamples ${TRAPDOORAPPS} )
	add_dependencies( alltrapdoor alltrapdoorexamples)
endif()

set( TRAPDOOREXTRAS "" )
if (BUILD_EXTRAS)
	file (GLOB TRAPDOOR_EXTRAS_SRC_FILES CONFIGURE_DEPENDS extras/*.cpp)
	foreach (app ${TRAPDOOR_EXTRAS_SRC_FILES})
		get_filename_component ( exe ${app} NAME_WE )
		add_executable ( ${exe} ${app} )
		set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/extras/trapdoor)
		set( TRAPDOOREXTRAS ${TRAPDOOREXTRAS} ${exe} )
		target_link_libraries ( ${exe} ${TRAPDOORLIBS} )
	endforeach()

	add_custom_target( alltrapdoorextras )
	add_dependencies( alltrapdoorextras ${TRAPDOOREXTRAS} )
	add_dependencies( alltrapdoor alltrapdoorextras)
endif()

add_custom_command( OUTPUT trapdoorinfocmd COMMAND echo Builds PALISADEtrapdoor and these apps: ${TRAPDOORAPPS} )
add_custom_target( trapdoorinfo DEPENDS trapdoorinfocmd )
