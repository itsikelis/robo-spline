# Copyright (c) 2023, Ioannis V. Tsikelis, University of Patras. All rights reserved.
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(namespace "rspl")
project("robo_spline" VERSION 1.0
    DESCRIPTION "A Cubic Spline Library fot Robot Trajectory Generation."
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include Eigen.
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Create executable
add_executable(${PROJECT_NAME} main.cpp)

# Link libraries and add compiler options
target_link_libraries(${PROJECT_NAME} PUBLIC
  Eigen3::Eigen
)

target_compile_options(${PROJECT_NAME} PUBLIC "-Wall;-Wextra;-O3") # -Werror

add_subdirectory(rspl)