cmake_minimum_required(VERSION 3.6)
project(FastTree)

set(CMAKE_CXX_STANDARD 11)

option(USE_NATIVE "enable/disable system's processor architecture optimization (linux)" OFF)
option(USE_SEE2 "enable/disable SEE2 in windows (linux default)" ON)
option(USE_SEE4 "enable/disable SEE4.1" OFF)
option(USE_AVX "enable/disable AVX" OFF)
option(USE_AVX2 "enable/disable AVX2" OFF)
option(USE_AVX512 "enable/disable AVX512" OFF)


add_definitions(-DNDEBUG)

ADD_SUBDIRECTORY("libs/CLI11")
ADD_SUBDIRECTORY("libs/robin-map")
ADD_SUBDIRECTORY("libs/xxhash")
ADD_SUBDIRECTORY("libs/boost-core")
ADD_SUBDIRECTORY("libs/boost-sort")
ADD_SUBDIRECTORY("libs/boost-align")

set(SOURCE_FILES
        fasttree/operations/BasicOperations.tcc
        fasttree/operations/BasicOperations.h
        fasttree/operations/SSE128Operations.tcc
        fasttree/operations/SSE128Operations.h
        fasttree/operations/AVX256Operations.tcc
        fasttree/operations/AVX256Operations.h
        fasttree/operations/AVX512Operations.tcc
        fasttree/operations/AVX512Operations.h
        fasttree/Knuth.cpp
        fasttree/Knuth.h
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

if (MSVC)
    set(CMAKE_CXX_FLAGS "/W2 /O2 /openmp /EHsc /bigobj")
    if (USE_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 /D__AVX2__=1 /D__AVX__=1")
    elseif (USE_AVX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX /D__AVX__=1")
    elseif (USE_SEE4)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2 /__SSE4_1__=1 /D__SSE2__=1")
    elseif (USE_SEE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2 /D__SSE2__=1")
    endif ()
else ()
    set(CMAKE_CXX_FLAGS "-Wall -fopenmp -O3")
    if (USE_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    elseif (USE_AVX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    elseif (USE_SEE4)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
    elseif (USE_AVX512)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512cd")
    endif ()
    if (USE_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif ()
endif ()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} CLI11 robin-map xxhash boost-core boost-sort boost-align)

