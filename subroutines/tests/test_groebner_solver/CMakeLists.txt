add_executable(test_groebner main.cpp simple_tests.cpp)

target_link_libraries(test_groebner PUBLIC test_utility PUBLIC libgtest PUBLIC libgmock)
target_include_directories(test_groebner PUBLIC ${PROJECT_SOURCE_DIR}/subroutines/distortion_groebner_estimator)
enable_testing()
#add_test(name testLinkedViewEstimator command testLVE "Test efficiency of rotation estimator")
