cmake_minimum_required(VERSION 3.16)
project(mecha)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")


include_directories(${PROJECT_SOURCE_DIR}/include/)
include_directories("D:/ftw/SENIOR/GNSSCenter/M5/LIB/eigen-3.3.9/")

add_subdirectory(${PROJECT_SOURCE_DIR}/src)


add_executable(mech_test mech_test.cpp)
target_link_libraries(mech_test Test)

add_executable(inter_test intergrated_test.cpp)
target_link_libraries(inter_test Test)