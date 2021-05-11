cmake_minimum_required(VERSION 3.15)
project(monte_carlo)

set(CMAKE_CXX_STANDARD 14)

add_executable(monte_carlo
        main.c
        makefile montecarlo.c montecarlo.h)
