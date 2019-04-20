project(robin-map)

add_library(robin-map INTERFACE)

target_include_directories(robin-map INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

