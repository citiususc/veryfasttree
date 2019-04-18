cmake_minimum_required(VERSION 3.6)
project(FastTree)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        libs/CLI11/CLI11.hpp
        libs/CLI11/Validators.h
        fasttree/Options.h
        fasttree/Constants.h
        fasttree/DebugVars.h
        fasttree/FastTree.cpp
        fasttree/FastTree.h
        fasttree/Utils.h
        main.cpp
        )

add_executable(FastTree ${SOURCE_FILES})
