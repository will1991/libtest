cmake_minimum_required(VERSION 2.8)

project(flann_test)
#find_package(LZ4 REQUIRED liblz4)
add_executable(${PROJECT_NAME} "main.cpp")
target_link_libraries(${PROJECT_NAME} lz4)
