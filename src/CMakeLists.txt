cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(cubic-hermite-splines VERSION 1.0)

set(CMAKE_DEBUG_POSTFIX d)

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_11)

set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
target_compile_options(compiler_flags INTERFACE
  "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused;-O3>>"
)

add_executable(cubic-hermite-splines main.cpp)

target_link_libraries(cubic-hermite-splines PUBLIC compiler_flags)
