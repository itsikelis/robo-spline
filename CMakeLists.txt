cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(cubic-hermite-splines VERSION 1.0)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Find all cpp files
file(GLOB_RECURSE CSPL_SRCS src/cspl/*.cpp)

# Create library
add_library(cspl ${CSPL_SRCS})
target_link_libraries(cspl PUBLIC Eigen3::Eigen)
target_include_directories(cspl PUBLIC src)
target_compile_options(cspl PUBLIC "-Wall;-Wextra;-Werror;-O3")

# Create example
add_executable(${PROJECT_NAME} src/examples/main.cpp)
# target_include_directories(${PROJECT_NAME} PUBLIC src)
target_link_libraries(${PROJECT_NAME} PUBLIC cspl)
target_compile_options(${PROJECT_NAME} PUBLIC "-Wall;-Wextra;-Werror;-O3")
