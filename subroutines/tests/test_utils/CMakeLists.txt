add_library(test_utility test_utils.h test_utils.cpp)

target_include_directories(test_utility PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})
