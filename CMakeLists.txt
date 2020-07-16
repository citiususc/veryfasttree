cmake_minimum_required(VERSION 3.5)
project(VeryFastTree)

set(CMAKE_CXX_STANDARD 11)

option(USE_NATIVE "enable/disable system's processor architecture optimization (linux)" ON)
option(USE_SEE2 "enable/disable SSE2 in windows (linux default)" ON)
option(USE_SEE4 "enable/disable SSE4.1" OFF)
option(USE_AVX "enable/disable AVX" OFF)
option(USE_AVX2 "enable/disable AVX2" OFF)
option(USE_AVX512 "enable/disable AVX512" OFF)


#add_definitions(-DNDEBUG)

ADD_SUBDIRECTORY("libs/CLI11")
ADD_SUBDIRECTORY("libs/robin-map")
ADD_SUBDIRECTORY("libs/xxhash")
ADD_SUBDIRECTORY("libs/boost-core")
ADD_SUBDIRECTORY("libs/boost-sort")
ADD_SUBDIRECTORY("libs/boost-align")

set(SOURCE_FILES
        src/operations/BasicOperations.tcc
        src/operations/BasicOperations.h
        src/operations/SSE128Operations.tcc
        src/operations/SSE128Operations.h
        src/operations/AVX256Operations.tcc
        src/operations/AVX256Operations.h
        src/operations/AVX512Operations.tcc
        src/operations/AVX512Operations.h
        src/Knuth.cpp
        src/Knuth.h
        src/Options.h
        src/Constants.h
        src/Debug.h
        src/Utils.h
        src/HashTable.h
        src/DistanceMatrix.tcc
        src/DistanceMatrix.h
        src/TransitionMatrix.tcc
        src/TransitionMatrix.h
        src/NeighbourJoining.tcc
        src/NeighbourJoining.h
        src/Alignment.cpp
        src/Alignment.h
        src/VeryFastTreeImpl.tcc
        src/VeryFastTree.cpp
        src/VeryFastTree.h
        main.cpp
        )

if (MSVC)
    set(CMAKE_CXX_FLAGS "/W2 /O2 /openmp /EHsc /bigobj")
    if (USE_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 /D__AVX2__=1 /D__AVX__=1")
    elseif (USE_AVX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX /D__AVX__=1")
    elseif (USE_SEE4)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2 /D__SSE4_1__=1 /D__SSE2__=1")
    elseif (USE_SEE2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2 /D__SSE2__=1")
    endif ()
else ()
    set(CMAKE_CXX_FLAGS "-Wall -fopenmp -O3")
    if (USE_AVX512)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f")
    elseif (USE_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
    elseif (USE_AVX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
    elseif (USE_SEE4)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
    endif ()
    if (USE_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif ()
endif ()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} CLI11 robin-map xxhash boost-core boost-sort boost-align)

