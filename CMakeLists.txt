cmake_minimum_required(VERSION 3.6)
project(FastTree)

set(CMAKE_CXX_STANDARD 11)

ADD_SUBDIRECTORY("libs/CLI11")
ADD_SUBDIRECTORY("libs/robin-map")

set(SOURCE_FILES
        fasttree/operations/BasicOperations.cpp
        fasttree/operations/BasicOperations.h
        fasttree/operations/SEEOperations.cpp
        fasttree/operations/SEEOperations.h
        fasttree/Options.h
        fasttree/Constants.h
        fasttree/Debug.h
        fasttree/Utils.h
        fasttree/FastTreeImpl.tcc
        fasttree/FastTreeImpl.h
        fasttree/FastTree.cpp
        fasttree/FastTree.h
        main.cpp
         )

if(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /O2")
else()
    set(CMAKE_CXX_FLAGS "-Wall")#-O3
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} CLI11 robin-map)

