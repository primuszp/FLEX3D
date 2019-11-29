cmake_minimum_required(VERSION 2.6)

# This is the 1-level cmake file for the whole project. Specify every directory
# of source code here.

# Project name
project(FEM)

# Specify directory(s) of source code
add_subdirectory(FEM)

# CMake tutorial
# http://www.hahack.com/codes/cmake/
# https://github.com/Campanula/CMake-tutorial

# Usage
# cd to the folder where the 1-level CMakeLists.txt (this file) locates. Under
# this folder, there should be a "/src" folder and a "/build" folder. "/src"
# contains all source code files, "/build" contains all compile-time files and
# the final target executable. If "/build" doesn't exist, "mkdir build" first.

# $ cd ./build && cmake .. // go into "/build" and build from the 1-level cmake config. This is called "out-of-source" build.
# $ make                   // Makefile will be generated under "/build". Now you can `make`
# $ ./src/main             // executable `main` is located under "/build/src/"
