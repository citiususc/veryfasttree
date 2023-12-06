cmake_minimum_required(VERSION 3.5...3.27)
project(bxzstr)


set(BXZSTR_Z_SUPPORT 1)
set(BXZSTR_BZ2_SUPPORT 1)
set(BXZSTR_LZMA_SUPPORT 0)
set(BXZSTR_ZSTD_SUPPORT 0)


configure_file(include/config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/config.hpp @ONLY)

add_library(bxzstr INTERFACE)
target_include_directories(bxzstr INTERFACE include)
