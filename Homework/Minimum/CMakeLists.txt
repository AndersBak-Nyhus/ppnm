cmake_minimum_required(VERSION 3.15)
project(minimization)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(minimization
        CMakeLists.txt
        main.c
        makefile
        minimization.c
        minimization.h)

find_package(GSL REQUIRED)    # See below (2)
target_link_libraries(minimization GSL::gsl GSL::gslcblas)