cmake_minimum_required(VERSION 3.5...3.27)
project(boost-align)

add_library(boost-align INTERFACE)

target_include_directories(boost-align INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

