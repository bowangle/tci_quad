cmake_minimum_required(VERSION 3.5)

project(TCIcpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ============================================================
# Include directory (for all targets)
# ============================================================
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/extern/eigen)

# ============================================================
# 1. Shared library (float128 instantiations)
# ============================================================
add_library(tci SHARED
    src/tci_quad_instantiations_float128.cpp
)

target_link_libraries(tci PRIVATE quadmath)

# ============================================================
# 2. Executable
# ============================================================
add_executable(TCIcpp
    src/main.cpp
)
add_executable(TestGrid
    test/test_grid.cpp
)

target_link_libraries(TestGrid PRIVATE quadmath)

target_link_libraries(TCIcpp PRIVATE quadmath tci)