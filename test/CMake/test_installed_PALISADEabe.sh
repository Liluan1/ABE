#!/bin/bash

function exit_with_msg() {
    echo "*********************************************************************"
    echo "${1}"
    echo "*********************************************************************"
    exit -1
}

# Build and install the PALISADEabe library
function build_lib() {
    rm -rf ../../build
    mkdir ../../build
    pushd ../../build 
    cmake .. || exit_with_msg "Failed to cmake config for lib"
    make -j16 || exit_with_msg "Failed to build lib"

    systemName=$(uname -s)
    case "${systemName}" in 
        Linux*)
            sudo make install || exit_with_msg "Failed to install lib"
            ;;
        MINGW*)     
            make install || exit_with_msg "Failed to install lib"
            ;;
        *)  
            exit_with_msg "system UNKNOWN:${systemName}"
            ;;
    esac
    popd
}
build_lib


# Build the examples with the installed PALISADEabe
function build_examples() {
    rm -rf build
    mkdir -p build
    pushd build
    cmake .. || exit_with_msg "Failed to config test harness"
    make -j16 || exit_with_msg "Failed to build examples with installed PALISADEabe package"
    popd
}
build_examples

# Run all build examples
for app in $(ls build/bin/); do
    ./build/bin/$app || exit "Failed ${app}"
done

echo "*********************************************************************"
echo "* Test Passed!"
echo "*********************************************************************"
