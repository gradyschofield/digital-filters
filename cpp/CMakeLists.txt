cmake_minimum_required(VERSION 3.19)
project(cpp)

set(CMAKE_CXX_STANDARD 20)
include_directories(src)

add_executable(cpp src/FilterOptimizer.cpp src/Polynomial.h src/Util.h src/DiskRational.h src/VecUtil.h)