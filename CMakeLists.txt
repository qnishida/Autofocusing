cmake_minimum_required(VERSION 3.16)
project(DeltaP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(DELTAP_ENABLE_NATIVE_ARCH "Build with -march=native on compilers that support it" OFF)

# コンパイラごとにオプションを設定
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Using GCC")
    set(COMPILER_SUFFIX "gcc")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "Using Clang")
    set(COMPILER_SUFFIX "clang")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(STATUS "Using MSVC")
    set(COMPILER_SUFFIX "msvc")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    message(STATUS "Using NVIDIA HPC (nvc++)")
    set(COMPILER_SUFFIX "nvc++")
else()
    message(WARNING "Unknown compiler: ${CMAKE_CXX_COMPILER_ID}")
    set(COMPILER_SUFFIX "unknown")
endif()

find_package(Boost REQUIRED COMPONENTS system filesystem date_time math_c99)
find_package(Eigen3 REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(OpenMP REQUIRED)

find_path(FFTW3_INCLUDE_DIR fftw3.h)
find_library(FFTW3_LIBRARY NAMES fftw3)
if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
    message(FATAL_ERROR "FFTW3 was not found. Set FFTW3_INCLUDE_DIR and FFTW3_LIBRARY if it is installed in a non-standard location.")
endif()

find_path(GEOGRAPHICLIB_INCLUDE_DIR GeographicLib/Geodesic.hpp)
find_library(GEOGRAPHICLIB_LIBRARY NAMES GeographicLib)
if(NOT GEOGRAPHICLIB_INCLUDE_DIR OR NOT GEOGRAPHICLIB_LIBRARY)
    message(FATAL_ERROR "GeographicLib was not found. Set GEOGRAPHICLIB_INCLUDE_DIR and GEOGRAPHICLIB_LIBRARY if it is installed in a non-standard location.")
endif()

set(DELTAP_INCLUDE_DIRS
    ${HDF5_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIR}
    ${GEOGRAPHICLIB_INCLUDE_DIR}
)

set(DELTAP_LIBS
    Boost::system
    Boost::filesystem
    Boost::date_time
    Boost::math_c99
    Eigen3::Eigen
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${FFTW3_LIBRARY}
    ${GEOGRAPHICLIB_LIBRARY}
    OpenMP::OpenMP_CXX
)

if(UNIX)
    list(APPEND DELTAP_LIBS m)
endif()

set(DELTAP_COMPILE_OPTIONS)
if(MSVC)
    list(APPEND DELTAP_COMPILE_OPTIONS /W4)
else()
    list(APPEND DELTAP_COMPILE_OPTIONS -Wall -Wextra -Wpedantic)
    if(DELTAP_ENABLE_NATIVE_ARCH AND
       (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
        list(APPEND DELTAP_COMPILE_OPTIONS -march=native)
    endif()
endif()

# テスト directory
set(TEST_OUTPUT_DIR ${CMAKE_BINARY_DIR}/tests)

# 実行ファイル名にコンパイラ名を付与
set(TARGET_NAME1 "cal_ccf_${COMPILER_SUFFIX}")
set(TARGET_NAME2 "cal_ccf_eq_${COMPILER_SUFFIX}")

add_subdirectory(src)

install(TARGETS ${TARGET_NAME1} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
