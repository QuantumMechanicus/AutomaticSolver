cmake_minimum_required(VERSION 3.8)
project(AutomaticSolver)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fsanitize=address,leak)
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,leak)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(subroutines)