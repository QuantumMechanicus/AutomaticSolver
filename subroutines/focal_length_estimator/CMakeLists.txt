FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
set(CMAKE_BUILD_TYPE Release)
set(SOURCE_FILES main_focal_length_estimator.cpp focal_length_solver.cpp focal_length_solver.h)
add_executable(FocalEstimator ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(FocalEstimator ${Boost_LIBRARIES})