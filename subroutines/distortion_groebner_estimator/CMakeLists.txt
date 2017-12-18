FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
set(CMAKE_BUILD_TYPE Release)
set(SOURCE_FILES main_distortion_estimator.cpp solver_ku8pt.cpp solver_ku8pt.h)
add_executable(GroebnerAutomaticSolver ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(GroebnerAutomaticSolver ${Boost_LIBRARIES} UndistortionUtils)