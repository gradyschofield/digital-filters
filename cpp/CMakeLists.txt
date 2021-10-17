cmake_minimum_required(VERSION 3.19)
project(cpp)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 23)


find_package(BLAS)

include_directories(src ${BLAS_INCLUDE_DIRS})

add_executable(cpp src/FilterOptimizer.cpp src/Polynomial.h src/Util.h src/DiskRational.h src/VecUtil.h src/BSpline.h src/Complex.h src/BLAS.h src/BLAS.cpp)
target_link_libraries(cpp ${BLAS_LIBRARIES})
