cmake_minimum_required(VERSION 3.5)

project(TCIcpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ============================================================
# Include directory (for all targets)
# ============================================================
set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(EIGEN_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/extern/eigen)
set(JSON_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/extern/nlohmann_json)

# ============================================================
# 1. Shared library (float128 instantiations)
# ============================================================
add_library(tci SHARED
    src/tci_quad_instantiations_float128.cpp
)

target_include_directories(tci PUBLIC
    ${PROJECT_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
    ${JSON_INCLUDE_DIR}
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
add_executable(TestTCI
    test/test_tci.cpp
)

target_include_directories(TCIcpp PRIVATE
    ${PROJECT_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
    ${JSON_INCLUDE_DIR}
)

target_include_directories(TestGrid PRIVATE
    ${PROJECT_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
    ${JSON_INCLUDE_DIR}
)

target_include_directories(TestTCI PRIVATE
    ${PROJECT_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
    ${JSON_INCLUDE_DIR}
)

target_link_libraries(TestGrid PRIVATE quadmath)
target_link_libraries(TestTCI PRIVATE quadmath)

target_link_libraries(TCIcpp PRIVATE quadmath tci)