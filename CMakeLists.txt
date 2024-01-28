cmake_minimum_required(VERSION 3.27)
project(PSO)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(PSO
        example/FindGlobalMinimum.cpp
        include/PSO.h
        include/BenchmarkFunction.h
)
