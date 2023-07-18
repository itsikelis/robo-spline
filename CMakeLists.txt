# Copyright (c) 2023, Ioannis V. Tsikelis, University of Patras. All rights reserved.
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(
    robo_spline
    VERSION 1.0
    DESCRIPTION "A Cubic Hermite Spline Library for Robot Trajectory Optimisation."
    LANGUAGES CXX)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    # Ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Support folders in IDEs
    # set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the
    # main CMakeLists.
    # include(CTest)

    # Docs only available if this is the main app
    find_package(Doxygen)
    if(Doxygen_FOUND)
      add_subdirectory(docs)
    else()
      message(STATUS "Doxygen not found, not building docs")
    endif()
endif()

# Include Eigen.
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(HEADER_LIST "${PROJECT_SOURCE_DIR}/include/robo_spline/trajectory.hpp")
set(INCLUDE_DEST "include/${PROJECT_NAME}${PROJECT_VERSION}")
set(LIB_DEST "lib/cmake/${PROJECT_NAME}${PROJECT_VERSION}")

add_subdirectory(src)

add_subdirectory(apps)

install(EXPORT ${PROJECT_NAME} DESTINATION "${LIB_DEST}")