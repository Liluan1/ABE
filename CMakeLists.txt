cmake_minimum_required (VERSION 3.5.1)

# this version will get incremented only if new changes are made in this repo
set(PALISADE_VERSION_MAJOR 1)
set(PALISADE_VERSION_MINOR 11)
set(PALISADE_VERSION_PATCH 0)
set(PALISADE_VERSION ${PALISADE_VERSION_MAJOR}.${PALISADE_VERSION_MINOR}.${PALISADE_VERSION_PATCH})

project(EXAMPLES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Palisade 1.11.0)

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
