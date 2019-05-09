cmake_minimum_required(VERSION 3.6)
project(FastTree)

set(CMAKE_CXX_STANDARD 11)

option(USE_SSE3 "enable/disable SSE3" OFF)
option(USE_AVX "enable/disable AVX" OFF)
option(USE_AVX2 "enable/disable AVX2" OFF)
option(USE_AVX512 "enable/disable AVX512" OFF)

if(${USE_SSE3})
    add_definitions(-DUSE_SSE3)
endif()

if(${USE_AVX})
    add_definitions(-DUSE_AVX)
endif()

if(${USE_AVX2})
    add_definitions(-DUSE_AVX2)
endif()

if(${USE_AVX512})
    add_definitions(-DUSE_AVX512)
endif()

ADD_SUBDIRECTORY("libs/CLI11")
ADD_SUBDIRECTORY("libs/robin-map")
ADD_SUBDIRECTORY("libs/xxhash")
ADD_SUBDIRECTORY("libs/boost-sort")
ADD_SUBDIRECTORY("libs/boost-align")

set(SOURCE_FILES
        fasttree/operations/BasicOperations.tcc
        fasttree/operations/BasicOperations.h
        fasttree/operations/SSEOperations.tcc
        fasttree/operations/SSEOperations.h
        fasttree/operations/AVXOperations.tcc
        fasttree/operations/AVXOperations.h
        fasttree/operations/AVX2Operations.tcc
        fasttree/operations/AVX2Operations.h
        fasttree/operations/AVX512Operations.tcc
        fasttree/operations/AVX512Operations.h
        fasttree/Options.h
        fasttree/Constants.h
        fasttree/Debug.h
        fasttree/Utils.h
        fasttree/HashTable.h
        fasttree/DistanceMatrix.tcc
        fasttree/DistanceMatrix.h
        fasttree/TransitionMatrix.tcc
        fasttree/TransitionMatrix.h
        fasttree/NeighbourJoining.tcc
        fasttree/NeighbourJoining.h
        fasttree/Alignment.cpp
        fasttree/Alignment.h
        fasttree/FastTreeImpl.tcc
        fasttree/FastTree.cpp
        fasttree/FastTree.h
        main.cpp
        )

if(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /O2 /openmp /EHsc /wd4100")
else()
    set(CMAKE_CXX_FLAGS "-Wall -fopenmp")#-O3
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} CLI11 robin-map xxhash boost-sort boost-align)

