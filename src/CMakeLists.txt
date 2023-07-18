add_library(${PROJECT_NAME} STATIC rspl.cpp ${HEADER_LIST})

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> # for headers when building
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # for config_impl.hpp when building
    $<INSTALL_INTERFACE:${INCLUDE_DEST}> # for client in install mode
    $<INSTALL_INTERFACE:${LIB_DEST}> # for config_impl.hpp in install mode
)

target_link_libraries(${PROJECT_NAME} PUBLIC
Eigen3::Eigen)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

target_compile_options(${PROJECT_NAME} PUBLIC "-Wall;-Wextra;-O3") # -Werror

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME} DESTINATION "${LIB_DEST}")
install(FILES ${HEADER_LIST} EXPORT ${PROJECT_NAME} DESTINATION "${INCLUDE_DEST}")  