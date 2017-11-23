cmake_minimum_required(VERSION 3.8)
find_package(OpenCV REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)

set(SOURCE_FILES main.cpp)


link_directories(${OpenCV_LIBS_DIR})
include_directories(${OpenCV_INCLUDE_PATH})

add_executable(Undistortion ${SOURCE_FILES})
target_link_libraries(Undistortion tbb ${OpenCV_LIBS} ${Boost_LIBRARIES})

