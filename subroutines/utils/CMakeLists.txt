add_library(UndistortionUtils STATIC undistortion_problem_utils.cpp undistortion_problem_utils.h)
target_include_directories(UndistortionUtils PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )