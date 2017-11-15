cmake_minimum_required(VERSION 3.8)
project(AutomaticSolver)

set(CMAKE_CXX_STANDARD 11)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp solver_ku8pt.cpp solver_ku8pt.h)
add_executable(AutomaticSolver ${SOURCE_FILES})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/pipeline)

TARGET_LINK_LIBRARIES(AutomaticSolver ${Boost_LIBRARIES})