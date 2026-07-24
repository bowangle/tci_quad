cmake_minimum_required(VERSION 3.16)
project(tci_quad LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)

# ---- dependency: QTgrid-quad manages itself (and transitively numeric-type) ---
add_subdirectory(extern/QTgrid-quad)

# ---- own paths ----------------------------------------------------------------
set(PROJECT_INCLUDE_DIR     ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(spdlog_INCLUDE_DIR      ${CMAKE_CURRENT_SOURCE_DIR}/extern/spdlog/include)

# ---- interface target: carries our headers + transitive deps ------------------
add_library(tci_quad_deps INTERFACE)
target_include_directories(tci_quad_deps INTERFACE
    ${PROJECT_INCLUDE_DIR}
    ${spdlog_INCLUDE_DIR}
)
target_link_libraries(tci_quad_deps INTERFACE qtgrid_deps OpenMP::OpenMP_CXX)

# ---- test / exec targets ------------------------------------------------------
add_executable(TCIcpp  src/main.cpp)
add_executable(TestTCI test/test_tci.cpp)

foreach(_tgt TCIcpp TestTCI)
    target_link_libraries(${_tgt} PRIVATE tci_quad_deps)
    target_compile_options(${_tgt} PRIVATE -Wall -Wextra -Wreorder -Wshadow)
endforeach()

