cmake_minimum_required (VERSION 3.5.1)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

project(PALISADEabe CXX)

# this version will get incremented only if new changes are made in this repo
set(ABE_VERSION_MAJOR 1)
set(ABE_VERSION_MINOR 11)
set(ABE_VERSION_PATCH 0)
set(ABE_VERSION ${ABE_VERSION_MAJOR}.${ABE_VERSION_MINOR}.${ABE_VERSION_PATCH})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Palisade 1.11.0 REQUIRED)
find_package(cereal 1.3.2 REQUIRED)

set( CMAKE_CXX_FLAGS ${PALISADE_CXX_FLAGS} )
set( CMAKE_EXE_LINKER_FLAGS ${PALISADE_EXE_LINKER_FLAGS} )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${PALISADE_INCLUDE} )
include_directories( ${PALISADE_INCLUDE}/third-party/include )
include_directories( ${PALISADE_INCLUDE}/core )
include_directories( ${PALISADE_INCLUDE}/pke )
include_directories( ${cereal_INCLUDE_DIRS})

#--------------------------------------------------------------------
# Compiler logic
#--------------------------------------------------------------------
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 6.1
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.1)
        message(WARNING "GCC version should be at least 6.1.")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    # require at least clang 6
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6)
        message(WARNING "Clang version should be at least 6.")
    endif()
else()
	message(WARNING "You are using ${CMAKE_CXX_COMPILER_ID} version ${CMAKE_CXX_COMPILER_VERSION}, which is unsupported.")
endif()

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${LIBINSTALL}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${LIBINSTALL}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${LIBINSTALL}")
endif("${isSystemDir}" STREQUAL "-1")

# Compiler flags

# Added -Wno-parentheses for compatibility with g++
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set (IGNORE_WARNINGS "-Wno-parentheses")
    # we can use GNU built-in functions provided by GCC for debugging. ex: __builtin_LINE (), __builtin_FUNCTION (), __builtin_FILE ()
    add_definitions(-DBUILTIN_INFO_AVAILABLE)
    message (STATUS "BUILTIN_INFO_AVAILABLE is defined")
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	set(IGNORE_WARNINGS "${IGNORE_WARNINGS} -Wno-unused-private-field -Wno-shift-op-parentheses")
endif()

set(C_COMPILE_FLAGS "-Wall -Werror -O3 ")
set(CXX_COMPILE_FLAGS "-Wall -Werror -Wno-reorder -O3 ${NATIVE_OPT} ${IGNORE_WARNINGS}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_COMPILE_FLAGS}")

#--------------------------------------------------------------------
# Installation logic
#--------------------------------------------------------------------
set(INSTALL_LIB_DIR lib CACHE PATH
	"Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include/palisade-abe CACHE PATH
	"Installation directory for headers")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/PALISADEabe)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

foreach(p LIB INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

message("***** INSTALL IS AT ${CMAKE_INSTALL_PREFIX}; to change, run cmake with -DCMAKE_INSTALL_PREFIX=/your/path")
set (CMAKE_INSTALL_MESSAGE LAZY)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

#--------------------------------------------------------------------
# Machine-specific checks
#--------------------------------------------------------------------
# determine the architecture on a Linux/Unix/macOS/MinGW system
if(CMAKE_HOST_UNIX OR MINGW)
  EXECUTE_PROCESS( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
else()
  set(ARCHITECTURE "unknown")
endif()

if(ARCHITECTURE)
  message(STATUS "Architecture is " ${ARCHITECTURE})
endif()

# third party directories

include_directories( third-party/google-test/googletest third-party/google-test/googletest/include )

### add directories for other PALISADE modules

link_directories( ${PALISADE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
link_libraries( ${PALISADE_LIBRARIES} )
message(STATUS "PALISADE_LIBRARIES found: ${PALISADE_LIBRARIES}")

option( BUILD_UNITTESTS "Set to ON to build unit tests for the library" ON)
option( BUILD_EXAMPLES "Set to ON to build examples for the library" ON)
option( GIT_SUBMOD_AUTO "Submodules auto-update" ON)

message( STATUS "BUILD_UNITTESTS:  ${BUILD_UNITTESTS}")
message( STATUS "BUILD_EXAMPLES:   ${BUILD_EXAMPLES}")
message( STATUS "GIT_SUBMOD_AUTO:  ${GIT_SUBMOD_AUTO}")

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

add_subdirectory(src/abe)

# auto-update git submodules
find_package (Git REQUIRED)
if(GIT_SUBMOD_AUTO AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    message(STATUS "Submodule update")
    if(NOT GIT_SUBMODULE_SYNCED)
        # "git submodule sync --recursive" should run only once, when CMakeCache.txt doesn't exist'
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule sync --recursive
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMODULE_RESULT)
        if(NOT GIT_SUBMODULE_RESULT EQUAL "0")
            message(FATAL_ERROR "\"git submodule sync --recursive\" failed with ${GIT_SUBMODULE_RESULT}, please checkout submodules")
        endif()
    endif()

    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE GIT_SUBMODULE_RESULT)
    if(NOT GIT_SUBMODULE_RESULT EQUAL "0")
        if(NOT GIT_SUBMODULE_SYNCED)
            # print this message only if update has never happened
            message(FATAL_ERROR "\"git submodule update --init\" failed with ${GIT_SUBMODULE_RESULT}, please checkout submodules")
        else()
            message(SEND_ERROR "\"git submodule update --init\" failed with ${GIT_SUBMODULE_RESULT}, please checkout submodules or disable autoupdate with -DGIT_SUBMOD_AUTO=OFF")
        endif()
    endif()

    if(NOT GIT_SUBMODULE_SYNCED)
        set(GIT_SUBMODULE_SYNCED ON CACHE BOOL "" FORCE)
    endif()
endif()

#--------------------------------------------------------------------
# Create PALISADEabe package
#--------------------------------------------------------------------
set(ABE_SHARED_LIBS "PALISADEabe")
set(ABE_PACKAGE_LIBS ${ABE_SHARED_LIBS})

export(EXPORT PALISADEabeTargets FILE "${PROJECT_BINARY_DIR}/PALISADEabeTargets.cmake")

export(PACKAGE PALISADEabe)

# Create the PALISADEabeConfig.cmake and PALISADEabeConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(PALISADEabeConfig.cmake.in
  "${PROJECT_BINARY_DIR}/PALISADEabeConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${PALISADE_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(PALISADEabeConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/PALISADEabeConfig.cmake" @ONLY)
# ... for both
configure_file(PALISADEabeConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/PALISADEabeConfigVersion.cmake" @ONLY)

# Install the PALISADEabeConfig.cmake and PALISADEabeConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/PALISADEabeConfig.cmake"
  "${PROJECT_BINARY_DIR}/PALISADEabeConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT PALISADEabeTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)


add_custom_target( allmodules DEPENDS ${PALISADE_SHARED_LIBS} )
