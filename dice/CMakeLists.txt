cmake_minimum_required(VERSION 3.9)
project(dice)

include_directories(${CMAKE_SOURCE_DIR}/deps)
set(CMAKE_CXX_STANDARD 11)

add_executable(dice main.cpp)
target_link_libraries(dice pthread)
