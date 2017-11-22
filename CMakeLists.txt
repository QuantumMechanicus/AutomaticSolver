cmake_minimum_required(VERSION 3.8)
project(AutomaticSolver)


set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fsanitize=address,leak)
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,leak)

set(CMAKE_CXX_STANDARD 11)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

find_package(Ceres REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${CERES_INCLUDE_DIRS})
set(SOURCE_FILES main.cpp solver_ku8pt.cpp solver_ku8pt.h)
add_executable(AutomaticSolver ${SOURCE_FILES})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/pipeline)

TARGET_LINK_LIBRARIES(AutomaticSolver ${Boost_LIBRARIES} ${CERES_LIBRARIES})