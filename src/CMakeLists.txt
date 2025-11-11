add_library(XFStream STATIC
	Nurbs.cxx
	LoopBuilder.cxx
	RegionBuilder.cxx
)

target_include_directories(XFStream PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(XFStream PUBLIC cxx_std_17)
