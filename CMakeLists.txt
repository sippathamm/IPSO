cmake_minimum_required(VERSION 3.27)
project(IPSO)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

add_executable(IPSO
        include/BenchmarkFunction.h
        include/IPSO.h
        example/TensionCompressionSpringDesign.cpp
)
