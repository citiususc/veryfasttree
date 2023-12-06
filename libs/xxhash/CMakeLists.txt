cmake_minimum_required(VERSION 3.5...3.27)
project(xxhash)

add_library(xxhash INTERFACE)

target_include_directories(xxhash INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

