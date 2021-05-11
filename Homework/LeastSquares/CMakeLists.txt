cmake_minimum_required(VERSION 3.15)
project(clion_testing C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(least-squares main.c leastSquares.c leastSquares.h gramSchmidt.c gramSchmidt.h input2Array.c input2Array.h utilities.c utilities.h)

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(least-squares GSL::gsl GSL::gslcblas)