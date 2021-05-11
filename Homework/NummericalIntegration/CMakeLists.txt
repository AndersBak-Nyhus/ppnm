cmake_minimum_required(VERSION 3.15)
project(numerical_integration)

set(CMAKE_CXX_STANDARD 14)

add_executable(numerical_integration
        main.c
        makefile integration.c integration.h)

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(numerical_integration GSL::gsl GSL::gslcblas)
