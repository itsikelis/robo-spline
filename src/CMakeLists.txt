set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

###########
## Build ##
###########
add_library(${LIB_CORE}
    trajectory.hpp
    cubic_hermite_polynomial.hpp
    cubic_hermite_polynomial_acc.hpp   
)
target_include_directories(${LIB_CORE} 
  PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/test>
    $<INSTALL_INTERFACE:include>
    ${EIGEN3_INCLUDE_DIR}
)
# better specify features you need, not c++ version
target_compile_features(${LIB_CORE} 
  PUBLIC 
    cxx_range_for 
    cxx_auto_type
)

#############
## Install ##
#############
# store in export-set for find-script generation
install(
  TARGETS ${LIB_CORE}
  EXPORT  ${LIB_CORE}-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
# Copy header files to usr/local/include/robo_spline/*
install(
  DIRECTORY include/robo_spline/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/robo_spline
)

# Generate the find scripts
install(
  EXPORT ${LIB_CORE}-targets
  FILE ${LIB_CORE}-targets.cmake
  NAMESPACE rspl::
  DESTINATION ${config_package_location}
)