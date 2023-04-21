cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(KT)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)
find_package(Boost REQUIRED)

add_definitions(-Wall -pedantic -g)

set(SOURCE_DIR "../src")
set(INCLUDE_DIR "../include")

include_directories(${INCLUDE_DIR})

add_executable(KT
    KT.cpp
    ${SOURCE_DIR}/gsl_interface.cpp
    ${SOURCE_DIR}/cauchy.cpp
    ${SOURCE_DIR}/phase.cpp
    ${SOURCE_DIR}/asymptotic.cpp
    ${SOURCE_DIR}/path.cpp
    ${SOURCE_DIR}/array.cpp
    ${SOURCE_DIR}/matching.cpp
    ${SOURCE_DIR}/splined_omnes.cpp
    ${SOURCE_DIR}/splined_path.cpp
    ${SOURCE_DIR}/path_eta_pi_pi.cpp
    ${SOURCE_DIR}/angular_average.cpp
    ${SOURCE_DIR}/dispersion_integral.cpp
    ${SOURCE_DIR}/basis_function.cpp
    ${SOURCE_DIR}/iterative_solution.cpp
    )
target_link_libraries(KT PRIVATE GSL::gsl GSL::gslcblas Boost::boost)

enable_testing()

add_test(NAME KT
    COMMAND ./KT)





