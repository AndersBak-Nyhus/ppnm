cmake_minimum_required(VERSION 3.15)
project(clion_testing C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(ODE main.c RK-ODE.c RK-ODE.h functions.c functions.h)

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(ODE GSL::gsl GSL::gslcblas)