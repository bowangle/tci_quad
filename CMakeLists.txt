cmake_minimum_required(VERSION 3.5)

project(TCIcpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(TCIcpp main.cpp
    lazymat.h
    prrrlu.h
    pimat.h
    mindex.h
    tensortrain.h
    tensorcross.h)

target_link_libraries(TCIcpp PRIVATE quadmath)
