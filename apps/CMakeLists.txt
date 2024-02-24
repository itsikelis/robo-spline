# Copyright (c) 2023, Ioannis V. Tsikelis, University of Patras. All rights reserved.

add_executable(test_regular test_regular.cpp)
target_compile_features(test_regular PRIVATE cxx_std_17)
target_link_libraries(test_regular PRIVATE ${PROJECT_NAME})

add_executable(test_acc test_acc.cpp)
target_compile_features(test_acc PRIVATE cxx_std_17)
target_link_libraries(test_acc PRIVATE ${PROJECT_NAME})