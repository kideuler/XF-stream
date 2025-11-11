cmake_minimum_required(VERSION 3.20)

project(XFStream VERSION 0.0.0 DESCRIPTION "XFStream library" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(src)

include(CTest)
enable_testing()

add_subdirectory(test)