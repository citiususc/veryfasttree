cmake_minimum_required(VERSION 3.5)
project(VeryFastTree)

set(CMAKE_CXX_STANDARD 11)

option(USE_SHARED "enable/disable linking of system shared libraries" OFF)
option(USE_NATIVE "enable/disable system's processor architecture optimization (linux)" ON)
option(USE_SEE2 "enable/disable SSE2 in windows (linux default)" ON)
option(USE_SEE4 "enable/disable SSE4.1" OFF)
option(USE_AVX "enable/disable AVX" OFF)
option(USE_AVX2 "enable/disable AVX2" OFF)
option(USE_AVX512 "enable/disable AVX512" OFF)
option(USE_CUDA "enable/disable CUDA" OFF)
set(CUDA_ARCH "80" CACHE STRING "change CUDA Architecture")

add_definitions(-DCLI11_BOOST_OPTIONAL=0)
add_definitions(-DNDEBUG)

ADD_SUBDIRECTORY("libs/CLI11")
ADD_SUBDIRECTORY("libs/bxzstr")
if (USE_SHARED AND NOT MSVC)
    set(LIBRARIES z bz2)
else ()
    ADD_SUBDIRECTORY("libs/boost-align")
    ADD_SUBDIRECTORY("libs/boost-core")
    ADD_SUBDIRECTORY("libs/boost-sort")
    ADD_SUBDIRECTORY("libs/bzip2")
    ADD_SUBDIRECTORY("libs/zlib")
    ADD_SUBDIRECTORY("libs/robin-map")
    ADD_SUBDIRECTORY("libs/xxhash")
    set(LIBRARIES boost-core boost-sort boost-align robin_map xxhash zlibstatic bzip2static)
endif ()

set(LIBRARIES CLI11 bxzstr ${LIBRARIES} )

set(SOURCE_FILES
        src/operations/BasicOperations.tcc
        src/operations/BasicOperations.h
        src/operations/SSE128Operations.tcc
        src/operations/SSE128Operations.h
        src/operations/AVX256Operations.tcc
        src/operations/AVX256Operations.h
        src/operations/AVX512Operations.tcc
        src/operations/AVX512Operations.h
        src/impl/VeryFastTreeDouble.cpp
        src/impl/VeryFastTreeDoubleAVX256.cpp
        src/impl/VeryFastTreeDoubleAVX512.cpp
        src/impl/VeryFastTreeDoubleCuda.cpp
        src/impl/VeryFastTreeDoubleSSE128.cpp
        src/impl/VeryFastTreeFloat.cpp
        src/impl/VeryFastTreeFloatAVX256.cpp
        src/impl/VeryFastTreeFloatAVX512.cpp
        src/impl/VeryFastTreeFloatCuda.cpp
        src/impl/VeryFastTreeFloatSSE128.cpp
        src/Knuth.cpp
        src/Knuth.h
        src/Options.h
        src/Constants.h
        src/Debug.h
        src/DiskMemory.cpp
        src/DiskMemory.h
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

if (USE_CUDA)
    message("WARNING!! CUDA is only a experimental feature")
    set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH})
    set(USE_CUDA 1)
    enable_language(CUDA)
    find_package(CUDAToolkit)
    set(SOURCE_FILES
            ${SOURCE_FILES}
            src/operations/CudaOperations.h
            src/operations/CudaOperations.cu
            )

    set(LIBRARIES ${LIBRARIES} CUDA::cublas)
else ()
    set(USE_CUDA 0)
endif ()

if (MSVC)
    ADD_SUBDIRECTORY("libs/mman-win32")
    set(LIBRARIES ${LIBRARIES} mman-win32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W2 /O2 /EHsc /bigobj /DUSE_CUDA=${USE_CUDA}")
    if (USE_AVX512)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512 /D__AVX512F__=1 /D__AVX2__=1 /D__AVX__=1")
    elseif (USE_AVX2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 /D__AVX2__=1 /D__AVX__=1")
    elseif (USE_AVX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX /D__AVX__=1")
    elseif (USE_SEE4)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2 /D__SSE4_1__=1 /D__SSE2__=1")
    elseif (USE_SEE2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2 /D__SSE2__=1")
    endif ()
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -DUSE_CUDA=${USE_CUDA}")
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

find_package(OpenMP)
if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX ${LIBRARIES})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)