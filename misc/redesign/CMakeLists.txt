# Copyright (c) 2023, Ioannis V. Tsikelis, University of Patras. All rights reserved.

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Declare Project Name, Descritption and Version.
project(
    robo_spline
    VERSION 1.0
    LANGUAGES CXX)

# Include Eigen.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(test_regular test_regular.cpp)
# add_executable(main2 main2.cpp)
target_compile_features(test_regular PRIVATE cxx_std_17)
target_link_libraries(test_regular PUBLIC Eigen3::Eigen)

# Set compiler flags.
target_compile_features(test_regular PRIVATE cxx_std_17)
target_compile_options(test_regular PRIVATE "-Wall;-Wextra;-O3") #-Werror;

add_executable(test_phased test_phased.cpp)
# add_executable(main2 main2.cpp)
target_compile_features(test_phased PRIVATE cxx_std_17)
target_link_libraries(test_phased PUBLIC Eigen3::Eigen)

# Set compiler flags.
target_compile_features(test_phased PRIVATE cxx_std_17)
target_compile_options(test_phased PRIVATE "-Wall;-Wextra;-O3") #-Werror;