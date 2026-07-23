cmake_minimum_required(VERSION 3.5)
project(TCIcpp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)

set(QTGRID_DIR              ${CMAKE_SOURCE_DIR}/extern/QTgrid-quad)
set(NUMERIC_TYPE_DIR        ${QTGRID_DIR}/extern/numeric-type-quad)

set(PROJECT_INCLUDE_DIR     ${CMAKE_SOURCE_DIR}/include)
set(EIGEN_INCLUDE_DIR       ${NUMERIC_TYPE_DIR}/extern/eigen)
set(QD_INCLUDE_DIR          ${NUMERIC_TYPE_DIR}/extern/qd-install/include)
set(QD_LIBRARY_DIR          ${NUMERIC_TYPE_DIR}/extern/qd-install/lib)
set(JSON_INCLUDE_DIR        ${QTGRID_DIR}/extern/json/include)
set(spdlog_INCLUDE_DIR      ${CMAKE_SOURCE_DIR}/extern/spdlog/include)
set(TYPING_INCLUDE_DIR      ${NUMERIC_TYPE_DIR}/include)
set(QTGRID_INCLUDE_DIR      ${QTGRID_DIR}/include)

add_executable(TCIcpp           src/main.cpp)
add_executable(TestTCI         test/test_tci.cpp)

# ---- include dirs (same for all targets) ----------------------------------------
foreach(_tgt TCIcpp TestTCI)
    target_include_directories(${_tgt} PUBLIC
        ${PROJECT_INCLUDE_DIR}
        ${EIGEN_INCLUDE_DIR}
        ${QD_INCLUDE_DIR}
        ${TCI_QUAD_DIR}
        ${JSON_INCLUDE_DIR}
        ${spdlog_INCLUDE_DIR}
        ${TYPING_INCLUDE_DIR}
        ${QTGRID_INCLUDE_DIR}
    )
    target_link_directories(${_tgt} PRIVATE ${QD_LIBRARY_DIR})
    target_link_libraries(${_tgt} PRIVATE qd quadmath OpenMP::OpenMP_CXX)
    target_compile_options(${_tgt} PRIVATE -Wall -Wextra -Wreorder -Wshadow)
endforeach()
