cmake_minimum_required(VERSION 3.16)
project(DeltaP LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(DELTAP_ENABLE_NATIVE_ARCH "Build with -march=native on compilers that support it" OFF)
option(DELTAP_ENABLE_METAL "Build the optional Apple Metal slant-stack backend" ${APPLE})
add_library(slant_stack_backend STATIC src/slant_stack_backend.cpp)
if(DELTAP_ENABLE_METAL)
    if(NOT APPLE)
        message(FATAL_ERROR "DELTAP_ENABLE_METAL requires Apple Metal; use OFF on other platforms")
    endif()
    enable_language(OBJCXX)
    target_sources(slant_stack_backend PRIVATE src/slant_stack_metal.mm)
    target_compile_definitions(slant_stack_backend PRIVATE AUTOFOCUSING_HAVE_METAL=1)
    target_compile_options(slant_stack_backend PRIVATE $<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-arc>)
    set_target_properties(slant_stack_backend PROPERTIES OBJCXX_STANDARD 14 OBJCXX_STANDARD_REQUIRED ON)
    find_library(METAL_FRAMEWORK Metal)
    find_library(FOUNDATION_FRAMEWORK Foundation)
    if(NOT METAL_FRAMEWORK OR NOT FOUNDATION_FRAMEWORK)
        message(FATAL_ERROR "Metal and Foundation frameworks are required")
    endif()
    target_link_libraries(slant_stack_backend PUBLIC ${METAL_FRAMEWORK} ${FOUNDATION_FRAMEWORK})
endif()

# コンパイラごとにオプションを設定
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Using GCC")
    set(COMPILER_SUFFIX "gcc")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
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

find_package(Boost REQUIRED COMPONENTS filesystem date_time math_c99)
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
    slant_stack_backend
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

include(CTest)
if(BUILD_TESTING)
    add_executable(test_io tests/io.cpp src/slant_stack.cpp src/station_info.cpp src/util.cpp src/calTT.cpp)
    target_include_directories(test_io PRIVATE ${DELTAP_INCLUDE_DIRS})
    target_link_libraries(test_io PRIVATE ${DELTAP_LIBS})
    add_test(NAME io_equivalence COMMAND test_io ${CMAKE_CURRENT_BINARY_DIR}/io-fixture.h5)
    if(DELTAP_ENABLE_METAL)
        add_executable(test_metal_slant_stack tests/metal_slant_stack.cpp src/slant_stack.cpp)
        target_link_libraries(test_metal_slant_stack PRIVATE slant_stack_backend OpenMP::OpenMP_CXX)
        # GPU access is required; run this executable explicitly on a Metal host.
    endif()
    add_executable(test_horizontal tests/horizontal.cpp src/slant_stack.cpp src/station_info.cpp src/util.cpp src/calTT.cpp)
    target_include_directories(test_horizontal PRIVATE ${DELTAP_INCLUDE_DIRS})
    target_link_libraries(test_horizontal PRIVATE ${DELTAP_LIBS})
    add_test(NAME horizontal_components COMMAND test_horizontal ${CMAKE_CURRENT_BINARY_DIR}/horizontal-fixture.h5)
    add_executable(test_slant_stack tests/slant_stack.cpp src/slant_stack.cpp src/station_info.cpp src/util.cpp src/calTT.cpp)
    target_include_directories(test_slant_stack PRIVATE ${DELTAP_INCLUDE_DIRS})
    target_link_libraries(test_slant_stack PRIVATE ${DELTAP_LIBS})
    add_test(NAME slant_stack_equivalence COMMAND test_slant_stack)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME parent_launcher COMMAND ${Python3_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/tests/launcher.py ${PROJECT_SOURCE_DIR}
            $<TARGET_FILE:${TARGET_NAME1}> ${CMAKE_CURRENT_BINARY_DIR})
    endif()
endif()
