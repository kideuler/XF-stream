cmake_minimum_required(VERSION 3.20)

project(XFStream VERSION 0.0.0 DESCRIPTION "XFStream library" LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# Required Gmsh support
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Gmsh REQUIRED)
message(STATUS "Gmsh found: ${GMSH_LIBRARY}")

# MFEM support
find_package(MFEM REQUIRED)
if(MFEM_FOUND)
	find_path(UMFPACK_INCLUDE_DIR
		NAMES umfpack.h
		HINTS /opt/homebrew/include /usr/local/include /usr/include
		PATH_SUFFIXES suitesparse SuiteSparse
	)
	if(UMFPACK_INCLUDE_DIR)
		message(STATUS "Found UMFPACK headers at: ${UMFPACK_INCLUDE_DIR}")
	else()
		message(STATUS "UMFPACK headers not found; if MFEM requires them, set UMFPACK_INCLUDE_DIR to your SuiteSparse include path.")
	endif()
endif()

add_subdirectory(src)

include(CTest)
enable_testing()

add_subdirectory(test)

# Examples
add_executable(make_shapes examples/make_shapes.cxx)
target_link_libraries(make_shapes PRIVATE XFStream)
target_compile_features(make_shapes PRIVATE cxx_std_17)

add_executable(make_mesh examples/make_mesh.cxx)
target_link_libraries(make_mesh PRIVATE XFStream)
target_compile_features(make_mesh PRIVATE cxx_std_17)

# Convenience target: generate meshes for example shapes with h=0.05
add_custom_target(example
	DEPENDS make_shapes make_mesh
	COMMAND $<TARGET_FILE:make_shapes>
	COMMAND $<TARGET_FILE:make_mesh> upper_half_disk.nbt 0.05 upper_half_disk.msh
	COMMAND $<TARGET_FILE:make_mesh> snail_shell.nbt 0.05 snail_shell.msh
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

