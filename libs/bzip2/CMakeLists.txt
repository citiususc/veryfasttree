cmake_minimum_required(VERSION 3.5...3.27)
project(bzip2 C)

if (MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /O2")
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
endif ()

set(BZLIB2_SRCS
    blocksort.c 
	huffman.c 
	crctable.c 
	randtable.c 
	compress.c 
	decompress.c 
	bzlib.c 
	bzip2recover.c
)


add_library(bzip2static STATIC ${BZLIB2_SRCS})
target_include_directories(bzip2static PUBLIC .)

