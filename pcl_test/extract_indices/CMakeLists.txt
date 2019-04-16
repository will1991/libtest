cmake_minimum_required(VERSION 2.8)
project(extract_indices)
find_package(PCL 1.8 REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})

add_executable(extract_indices extract_indices.cpp )
target_link_libraries(extract_indices ${PCL_LIBRARIES})