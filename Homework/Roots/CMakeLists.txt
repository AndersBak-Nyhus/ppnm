cmake_minimum_required(VERSION 3.15)
project(root-finding)

set(CMAKE_CXX_STANDARD 14)

add_executable(root-finding
        main.c
        makefile rootfinding.c rootfinding.h gramSchmidt.h gramSchmidt.c RK-ODE.c RK-ODE.h)

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(root-finding GSL::gsl GSL::gslcblas)