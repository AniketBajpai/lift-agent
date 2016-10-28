cmake_minimum_required(VERSION 3.4)
project(lift_agent)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES State.cpp Problem.cpp)
add_executable(lift_agent ${SOURCE_FILES})