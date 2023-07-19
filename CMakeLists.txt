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

    # # Docs only available if this is the main app
    # find_package(Doxygen)
    # if(Doxygen_FOUND)
    #   add_subdirectory(docs)
    # else()
    #   message(STATUS "Doxygen not found, not building docs")
    # endif()
endif()

# Include Eigen.
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include(GNUInstallDirs)

set(config_package_location "share/${PROJECT_NAME}/cmake") # for .cmake find-scripts installs

add_subdirectory(src)

add_subdirectory(apps)

# generate the file IfoptConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
# install the two files in a place where cmake looks for them
install(
  FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}-config.cmake"         # self-written
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake" # generated
  DESTINATION ${config_package_location}
)

install(
  EXPORT ${PROJECT_NAME}-targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${config_package_location}
)

# uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

