# Copyright (c) 2023, Ioannis V. Tsikelis, University of Patras. All rights reserved.

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# Declare Project Name, Descritption and Version.
project(
    robo_spline
    VERSION 1.0
    DESCRIPTION "A Cubic Hermite Spline Library for Robot Trajectory Optimisation."
    LANGUAGES CXX)

# Generate Docs
# find_package(Doxygen)
# if(Doxygen_FOUND)
#   add_subdirectory(docs)
# else()
#   message(STATUS "Doxygen not found, not building docs")
# endif()

# Include Eigen.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Include GNU Install Directory Variables.
include(GNUInstallDirs)

# Location for CMake to search for Config file.
set(CONFIG_PACKAGE_LOCATION "share/${PROJECT_NAME}/cmake")

add_subdirectory(src)

add_subdirectory(apps)

# Generate the robo_spline-config-version.cmake file.
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# Install Config (Lists Library dependencies) and Config Version files.
install(
  FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake"        # Self-written.
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake" # Generated.
  DESTINATION ${CONFIG_PACKAGE_LOCATION}
)

# Install Targets file.
install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CONFIG_PACKAGE_LOCATION}
)

# Add an Uninstall Target to uninstall Library.
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Uninstall.cmake)
endif()

