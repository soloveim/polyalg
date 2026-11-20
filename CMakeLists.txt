cmake_minimum_required(VERSION 3.14)
project(polyalg VERSION 1.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(POLYALG_BUILD_EXAMPLES "Build example programs" ON)
option(POLYALG_BUILD_TESTS "Build unit tests" ON)

# Header-only library target (polynomial.hpp, quadrature.hpp)
add_library(polyalg INTERFACE)
add_library(polyalg::polyalg ALIAS polyalg)
target_include_directories(polyalg INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_features(polyalg INTERFACE cxx_std_17)

function(polyalg_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endfunction()

if(POLYALG_BUILD_EXAMPLES)
    add_executable(polyalg_demo src/main.cpp)
    add_executable(quadrature_demo src/quadrature_demo.cpp)
    foreach(t polyalg_demo quadrature_demo)
        target_link_libraries(${t} PRIVATE polyalg::polyalg)
        polyalg_warnings(${t})
    endforeach()
endif()

if(POLYALG_BUILD_TESTS)
    enable_testing()
    add_executable(polyalg_tests tests/test_polynomial.cpp)
    add_executable(quadrature_tests tests/test_quadrature.cpp)
    foreach(t polyalg_tests quadrature_tests)
        target_link_libraries(${t} PRIVATE polyalg::polyalg)
        polyalg_warnings(${t})
        add_test(NAME ${t} COMMAND ${t})
    endforeach()
endif()

install(DIRECTORY include/ DESTINATION include)
