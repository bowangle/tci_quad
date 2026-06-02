cmake_minimum_required(VERSION 3.5)

project(TCIcpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ============================================================
# 1. Shared library (float128 instantiations)
# ============================================================

add_library(tci SHARED
    tci_quad_instantiations_float128.cpp
)

target_include_directories(tci PRIVATE
    ${CMAKE_SOURCE_DIR}/extern/eigen
)

target_link_libraries(tci PRIVATE quadmath)

# ============================================================
# 2. Executable
# ============================================================

add_executable(TCIcpp
    main.cpp
    lazymat.h
    prrrlu.h
    pimat.h
    mindex.h
    tensortrain.h
    tensorcross.h
)

target_include_directories(TCIcpp PRIVATE
    ${CMAKE_SOURCE_DIR}/extern/eigen
)

target_link_libraries(TCIcpp PRIVATE quadmath tci)