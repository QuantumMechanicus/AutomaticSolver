FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

FIND_PACKAGE(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${CERES_INCLUDE_DIRS})
set(CMAKE_BUILD_TYPE Release)
set(SOURCE_FILES main_non_linear_optimizator.cpp)
add_executable(NonLinearOptimizator ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(NonLinearOptimizator ${Boost_LIBRARIES} ${CERES_LIBRARIES} undistortion_utils)