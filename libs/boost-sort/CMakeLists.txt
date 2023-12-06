cmake_minimum_required(VERSION 3.5...3.27)
project(boost-sort)

add_library(boost-sort INTERFACE)

target_include_directories(boost-sort INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

