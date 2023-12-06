cmake_minimum_required(VERSION 3.5...3.27)
project(mman-win32 C)

add_library(mman-win32 STATIC src/sys/mman.c)
target_include_directories(mman-win32 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
