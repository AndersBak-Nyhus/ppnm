cmake_minimum_required(VERSION 3.15)
project(clion_testing C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(eigenvalues main.c jacobi.c jacobi.h utilities.h utilities.c diffClock.c diffClock.h test_runtime.c)

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(eigenvalues GSL::gsl GSL::gslcblas)