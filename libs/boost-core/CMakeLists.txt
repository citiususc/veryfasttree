cmake_minimum_required(VERSION 3.5...3.27)
project(boost-core)

add_library(boost-core INTERFACE)

target_include_directories(boost-core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

