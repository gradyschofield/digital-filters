cmake_minimum_required(VERSION 3.19)
project(cpp)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 23)

find_package(BLAS)

include_directories(src ${BLAS_INCLUDE_DIRS})

set(SRC
        src/Polynomial.h
        src/Util.h
        src/DiskRational.h
        src/VecUtil.h
        src/BSpline.h
        src/Complex.h
        src/BLAS.h
        src/BLAS.cpp
        src/DerivativeOptimizable.h
        src/LineSearch.h
        src/GeneticAlgorithm.h
        src/GeneticallyOptimizable.h
        )

add_executable(filterOptimizer src/FilterOptimizer.cpp ${SRC} )
target_link_libraries(filterOptimizer ${BLAS_LIBRARIES})

add_executable(filterPeakFinder src/FilterPeakFinder.cpp ${SRC})
target_link_libraries(filterPeakFinder ${BLAS_LIBRARIES})
