add_library(undistortion_utils STATIC undistortion_problem_utils.cpp undistortion_problem_utils.h)
target_include_directories(undistortion_utils PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )