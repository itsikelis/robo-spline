# Copyright (c) 2023, Ioannis V. Tsikelis, University of Patras. All rights reserved.

# Include auto-generated Targets file.
include("${CMAKE_CURRENT_LIST_DIR}/robo_splineTargets.cmake")

get_target_property(robo_spline_INCLUDE_DIRS robo_spline::robo_spline INTERFACE_INCLUDE_DIRECTORIES)

get_property(robo_spline_LIB_CORE TARGET robo_spline::robo_spline PROPERTY LOCATION)
list(APPEND robo_spline_LIBRARIES ${robo_spline})
