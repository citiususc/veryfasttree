project(xxhash)

add_library(xxhash INTERFACE)

target_include_directories(xxhash INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

