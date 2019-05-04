project(CLI11)

add_library(boost-sort INTERFACE)

target_include_directories(boost-sort INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

