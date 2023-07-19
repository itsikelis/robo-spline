add_executable(rspl_test rspl_test.cpp)
target_compile_features(rspl_test PRIVATE cxx_std_17)

target_link_libraries(rspl_test PRIVATE ${PROJECT_NAME})
