add_library(${PROJECT_NAME} STATIC rspl.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> # for headers when building
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # for config_impl.hpp when building
)

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

target_compile_options(${PROJECT_NAME} PUBLIC "-Wall;-Wextra;-O3") # -Werror

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/include/robo_spline/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/robo_spline
)
