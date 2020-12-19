add_custom_target( allabe ) #special make target to make all abe files

include_directories(include) #to find abe/include files

#add_subDirectory(lib)


##########################################################
### build PALISADEabe share library and put in bin/lib
##########################################################


#list all sources in thie directory for the library (found in the lib subdirectory)
file(GLOB ABELIB_SOURCES "lib/*-impl.cpp")

#generate the shared library from the sources in that direcory
message(STATUS "{ABELIB_SOURCES} ${ABELIB_SOURCES}")
add_library(PALISADEabe SHARED ${ABELIB_SOURCES})

set_property(TARGET PALISADEabe PROPERTY POSITION_INDEPENDENT_CODE 1)

#set the local destination target for library output
#note without this we would build to build/src/lib
set_property(TARGET PALISADEabe PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

install(TARGETS PALISADEabe EXPORT PalisadeTargets DESTINATION lib)

#target_compile_options(pre-server-sockets PUBLIC -DUSE_SOCKETS)

add_dependencies( allabe PALISADEabe) #add library to make allabe target

##########################################################
### now build all files in examples and put executables in bin/abe/examples
##########################################################

#note PALISADEabe is our library build in lib

#set the libraries needed for our examples
set (ABELIBS PUBLIC PALISADEabe PUBLIC PALISADEcore ${THIRDPARTYLIBS} ${OpenMP_CXX_FLAGS})

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

#to add compiler options use: 
#target_compile_options(pre-server-sockets PUBLIC -DUSE_SOCKETS)


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


