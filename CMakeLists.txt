cmake_minimum_required(VERSION 3.6)
project(FastTree)

set(CMAKE_CXX_STANDARD 11)

ADD_SUBDIRECTORY("libs/CLI11")
ADD_SUBDIRECTORY("libs/robin-map")

set(SOURCE_FILES
        fasttree/Options.h
        fasttree/Constants.h
        fasttree/DebugVars.h
        fasttree/FastTree.cpp
        fasttree/FastTree.h
        fasttree/Utils.h
        main.cpp
        )

if(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /O2")
else()
    set(CMAKE_CXX_FLAGS "-Wall")#-O3
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} CLI11 robin-map)

