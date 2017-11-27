cmake_minimum_required(VERSION 3.8)
project(AutomaticSolver)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/pipeline)

#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fsanitize=address,leak)
#set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,leak)
set(CMAKE_CXX_STANDARD 14)


add_subdirectory(subroutines)
#add_subdirectory(subroutines/distortionEstimationGroebner)
#add_subdirectory(subroutines/nonLinearOptimization)
#add_subdirectory(subroutines/undistortion)

