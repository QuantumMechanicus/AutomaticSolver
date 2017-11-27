find_package(OpenCV REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)


link_directories(${OpenCV_LIBS_DIR})
include_directories(${OpenCV_INCLUDE_PATH})

add_executable(Undistortion main.cpp)
target_link_libraries(Undistortion tbb ${OpenCV_LIBS} ${Boost_LIBRARIES})

