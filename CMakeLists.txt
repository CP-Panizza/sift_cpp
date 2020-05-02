cmake_minimum_required(VERSION 3.12)
project(sift_cpp)

set(CMAKE_CXX_STANDARD 14)
link_libraries(C:\\Users\\Administrator\\Desktop\\sift_cpp\\img.a)
add_executable(sift_cpp
        gaussfilter.hpp
        img.h
        main.cpp
        matrix.hpp
        pgmer.hpp
        sift.hpp
        utils.cpp
        utils.h)
