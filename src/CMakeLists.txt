add_library(XFStream SHARED
	Nurbs.cxx
	LoopBuilder.cxx
	RegionBuilder.cxx
	MeshGenerator.cxx
	NbsIO.cxx
	NbtIO.cxx
	Mesh.cxx
	Laplace.cxx
)

target_include_directories(XFStream PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(XFStream PUBLIC cxx_std_17)
target_link_libraries(XFStream PUBLIC Gmsh::gmsh)

# Link MFEM if found (from root find_package)
if(MFEM_FOUND)
	# MFEM CMake export provides target 'MFEM::mfem' in recent versions; fall back to variables otherwise
	if(TARGET MFEM::mfem)
		target_link_libraries(XFStream PUBLIC MFEM::mfem)
	elseif(TARGET mfem)
		target_link_libraries(XFStream PUBLIC mfem)
	elseif(MFEM_LIBRARY)
		target_include_directories(XFStream PUBLIC ${MFEM_INCLUDE_DIRS})
		target_link_libraries(XFStream PUBLIC ${MFEM_LIBRARY})
	elseif(MFEM_LIBRARIES)
		target_include_directories(XFStream PUBLIC ${MFEM_INCLUDE_DIRS})
		target_link_libraries(XFStream PUBLIC ${MFEM_LIBRARIES})
	endif()
	# Homebrew SuiteSparse headers default path on Apple Silicon
	find_path(SUITESPARSE_INCLUDE_DIR
		NAMES umfpack.h cholmod.h SuiteSparse_config.h
		HINTS /opt/homebrew/include /usr/local/include /usr/include
		PATH_SUFFIXES suitesparse SuiteSparse
	)
	if(SUITESPARSE_INCLUDE_DIR)
		target_include_directories(XFStream PUBLIC ${SUITESPARSE_INCLUDE_DIR})
	endif()
endif()
