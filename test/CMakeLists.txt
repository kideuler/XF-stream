# Use system-installed GoogleTest
find_package(GTest REQUIRED)

# Support both package-config and FindGTest module target names
if (TARGET GTest::gtest_main)
	set(GTEST_MAIN_TARGET GTest::gtest_main)
elseif (TARGET GTest::Main)
	set(GTEST_MAIN_TARGET GTest::Main)
else()
	message(FATAL_ERROR "Could not find GTest main target (GTest::gtest_main or GTest::Main)")
endif()

if (TARGET GTest::gtest)
	set(GTEST_LIB_TARGET GTest::gtest)
elseif (TARGET GTest::GTest)
	set(GTEST_LIB_TARGET GTest::GTest)
else()
	set(GTEST_LIB_TARGET "")
endif()

add_executable(test_nurbs
	test_nurbs.cxx
	test_loop_builder.cxx
	test_region_builder.cxx
)
target_link_libraries(test_nurbs PRIVATE ${GTEST_MAIN_TARGET} ${GTEST_LIB_TARGET} XFStream)
target_compile_features(test_nurbs PRIVATE cxx_std_17)

include(GoogleTest)
gtest_discover_tests(test_nurbs)
