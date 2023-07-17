set(HEADER_LIST "${RoboSpline_SOURCE_DIR}/include/robo_spline/trajectory.hpp")

add_library(${PROJECT_NAME} STATIC rspl.cpp ${HEADER_LIST})

target_include_directories(${PROJECT_NAME} PUBLIC ../include)

target_link_libraries(${PROJECT_NAME} PUBLIC
Eigen3::Eigen)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

target_compile_options(${PROJECT_NAME} PUBLIC "-Wall;-Wextra;-O3") # -Werror