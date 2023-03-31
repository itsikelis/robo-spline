cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(cubic-hermite-splines VERSION 1.0)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(${PROJECT_NAME} src/main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

target_compile_options(${PROJECT_NAME} PUBLIC "-Wall;-O3")
