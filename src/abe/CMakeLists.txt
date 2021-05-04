#list all sources in this directory for the library (found in the lib subdirectory)
file(GLOB ABELIB_SOURCES "lib/*-impl.cpp")

include_directories(${CORE_INCLUDE_DIRS})
include_directories(${PKE_INCLUDE_DIRS})
list(APPEND ABE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include")
list(APPEND ABE_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib")
include_directories(${ABE_INCLUDE_DIRS})

add_library (abeobj OBJECT ${ABELIB_SOURCES})
add_dependencies(abeobj PALISADEcore)
set_property(TARGET abeobj PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(PALISADEabe SHARED $<TARGET_OBJECTS:abeobj>)
set_property(TARGET PALISADEabe PROPERTY VERSION ${ABE_VERSION})
set_property(TARGET PALISADEabe PROPERTY SOVERSION ${ABE_VERSION_MAJOR})
set_property(TARGET PALISADEabe PROPERTY POSITION_INDEPENDENT_CODE 1)

#sets the local destination target for library output
#this creates the lib directory in build and puts PALISADEabe there
#note without this we would build to build/src/lib
# this line is for non-Windows systems
set_property(TARGET PALISADEabe PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# this line is for Windows (MinGW) systems
set_property(TARGET PALISADEabe PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#these lines are needed for make install to work

install(
	TARGETS PALISADEabe 
	EXPORT PALISADEabeTargets 
	DESTINATION lib)
install(
	DIRECTORY include/
	DESTINATION include/palisade-abe)

set (ABELIBS PUBLIC PALISADEabe ${PALISADE_LIBRARIES})

target_link_libraries (PALISADEabe ${PALISADE_LIBRARIES})

add_custom_target( allabe ) #special make target to make all abe files
add_dependencies( allabe PALISADEabe) #add library to make allabe target

##########################################################
### now build all files in examples and put executables in bin/abe/examples
##########################################################

set (ABEAPPS "")
if( BUILD_EXAMPLES ) #set during cmake with -D
  #for all src files in examples directory
  file (GLOB ABE_EXAMPLES_SRC_FILES CONFIGURE_DEPENDS examples/*.cpp)
  foreach (app ${ABE_EXAMPLES_SRC_FILES})
	get_filename_component ( exe ${app} NAME_WE )
	add_executable ( ${exe} ${app} )
	##this is out output directory under build
	set_property(TARGET ${exe} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/examples/abe)
	set( ABEAPPS ${ABEAPPS} ${exe} ) ##append this apname
	target_link_libraries ( ${exe} ${ABELIBS} )
  endforeach()
  
  add_custom_target( allabeexamples ) #new make target
  add_dependencies( allabeexamples ${ABEAPPS} ) #builds all example apps
  add_dependencies( allabe allabeexamples ) 
endif()

##########################################################
### now build unittests
##########################################################

if( BUILD_UNITTESTS )
	set(UNITTESTMAIN ${PROJECT_SOURCE_DIR}/test/Main_TestAll.cpp)

  file (GLOB ABE_TEST_SRC_FILES CONFIGURE_DEPENDS unittest/*.cpp)
	add_executable (abe_tests ${ABE_TEST_SRC_FILES} ${UNITTESTMAIN})
	set_property(TARGET abe_tests PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/unittest)
	target_link_libraries ( abe_tests ${ABELIBS} )
	if (NOT ${USE_OpenMP} )
		target_link_libraries ( abe_tests PRIVATE Threads::Threads)
	endif()
	add_dependencies( allabe abe_tests )

	add_custom_command( OUTPUT runabetests WORKING_DIRECTORY ${CMAKE_BINARY_DIR} COMMAND ${CMAKE_BINARY_DIR}/unittest/abe_tests )
	add_custom_target( testabe DEPENDS abe_tests runabetests )
endif()

add_custom_command( OUTPUT abeinfocmd COMMAND echo Builds PALISADEabe and these apps: ${PKEAPPS} )
add_custom_target( abeinfo DEPENDS abeinfocmd )
