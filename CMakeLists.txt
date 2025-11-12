cmake_minimum_required(VERSION 3.20)

project(XFStream VERSION 0.0.0 DESCRIPTION "XFStream library" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# Required Gmsh support
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Gmsh REQUIRED)
message(STATUS "Gmsh found: ${GMSH_LIBRARY}")

add_subdirectory(src)

include(CTest)
enable_testing()

add_subdirectory(test)