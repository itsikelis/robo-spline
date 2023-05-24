# Copyright (c) 2023, Ioannis V. Tsikelis, University of Patras. All rights reserved.
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(robo_spline VERSION 1.0)

# find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# # Find all cpp files
# file(GLOB_RECURSE CSPL_SRCS src/cspl/*.cpp)

# # Create library
# add_library(cspl ${CSPL_SRCS})
# target_link_libraries(cspl PUBLIC Eigen3::Eigen)
# target_include_directories(cspl PUBLIC src)
# target_compile_options(cspl PUBLIC "-Wall;-Wextra;-Werror;-O3")

# # Create example
# add_executable(${PROJECT_NAME} src/examples/main.cpp)
# # target_include_directories(${PROJECT_NAME} PUBLIC src)
# target_link_libraries(${PROJECT_NAME} PUBLIC cspl)
# target_compile_options(${PROJECT_NAME} PUBLIC "-Wall;-Wextra;-Werror;-O3")


###########
## Build ##
###########
include(GNUInstallDirs) # for correct library locations across platforms
set(config_package_location "share/${PROJECT_NAME}/cmake") # for .cmake find-scripts installs
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)             # so installed solver libraries link to IFOPT/SNOPT                                          
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)    # so installed solver libraries link to libifopt_core.so  

IF(NOT CMAKE_BUILD_TYPE MATCHES Release)
  message(STATUS "CMAKE_BUILD_TYPE not set to Release -> impacts performance")
endif()

add_subdirectory(src)

##################################################
## Install find scripts for find_package(robo_spline) ##
##################################################
# generate the file IfoptConfigVersion.cmake
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/robo-spline-config-version.cmake"
  VERSION ${robo_spline_VERSION}
  COMPATIBILITY SameMajorVersion
)
# install the two files in a place where cmake looks for them
install(
  FILES 
    "${CMAKE_CURRENT_SOURCE_DIR}/robo-spline-config.cmake"         # self-written
    "${CMAKE_CURRENT_BINARY_DIR}/robo-spline-config-version.cmake" # generated
  DESTINATION ${config_package_location}
)

##########################
## Build binary package ##
##########################
if(BUILD_PACKAGE)
  set(DESCRIPTION "A Cubic Spline Library fot Robot Trajectory Generation.")
  set(LINUX_DEPENDS "libeigen3-dev")

  set(CPACK_PACKAGE_VENDOR "Ioannis V. Tsikelis")
#   set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_LIST_DIR}/LICENSE)
  set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_LIST_DIR}/README.md)
  if(UNIX)
    set(CPACK_GENERATOR "DEB;TXZ")

    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
      set(DEB_ARCH "amd64")
    else()
      set(DEB_ARCH ${CMAKE_SYSTEM_PROCESSOR})
    endif()

    set(CPACK_PACKAGE_FILE_NAME "${ARG_PACKAGE_PREFIX}${PROJECT_NAME}_${DEB_ARCH}_linux_${CMAKE_PROJECT_VERSION}")
    set(CPACK_DEBIAN_PACKAGE_NAME "${ARG_PACKAGE_PREFIX}${PROJECT_NAME}")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${DEB_ARCH})
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER <https://github.com/itsikelis/robo-spline>)
    set(CPACK_DEBIAN_PACKAGE_DESCRIPTION ${DESCRIPTION})
    set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS=ON)
    string(REPLACE ";" "," CPACK_DEBIAN_PACKAGE_DEPENDS "${LINUX_DEPENDS}")
  include(CPack)
  endif()
endif()