#!/bin/bash
# will manually build Catchment-CN according to the instructions from Matt
# Thompson

# check to see if the file exists
FILE=./build/
if [ -d "$FILE" ]; then
    # ask the user if we should make clean
    read -rep $'Make clean? <y/n>\n' clean
    # if make clean, remove build directory and recreate it
    if [ $clean = y ]; then
        echo "running make clean"
        rm -rf ./build
        mkdir ./build
        cd ./build
        cmake .. -DCMAKE_INSTALL_PREFIX=../install -DUSE_F2PY=OFF |& tee cmake.log
        make -j 6 install |& tee make.log
    # else do nothing
    else
        echo "not running make clean"
        cd ./build
        make install
    fi
# else if file does not exist then just create it and run make install from
# that directory
else
    echo "creating build directory and building model"
    mkdir ./build/
    cd ./build
    cmake .. -DCMAKE_INSTALL_PREFIX=../install -DUSE_F2PY=OFF |& tee cmake.log
    make -j 6 install |& tee make.log
fi
# exit out of build directory for all cases
cd ../
