#!/bin/bash
# bash script will run as in a sandbox, so cd here won't change the pwd outside
# remember to `chmod +x compile.sh` for this file or use `sh compile.sh`

# 1. setup build directory
[[ -d build ]] || mkdir build

# rm build/CMakeCache.txt # if the directory is moved or reconstructed, reset the cmake cache before make. When you see error, just remove this cache file

# 2. generate Makefile
cd ./build
cmake ..

# 3. make executable
make
# executable will be under ./build/FEM/main
