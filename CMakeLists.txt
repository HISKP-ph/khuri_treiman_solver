cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
cmake_policy(SET CMP0127 NEW)
project(khuri)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)

add_definitions(-Wall -Werror -pedantic -O3)

set(SOURCE_DIR "./cpp/src")
set(INCLUDE_DIR "./cpp/include")
set(BINDING_DIR "./khuri/bindings")

include_directories(${EIGEN3_INCLUDE_DIRS} ${INCLUDE_DIR})

add_subdirectory(lib/pybind11)

pybind11_add_module(_khuri_gsl
    "${SOURCE_DIR}/gsl_interface.cpp"
    "${BINDING_DIR}/gsl_bindings.cpp")
target_link_libraries(_khuri_gsl PRIVATE GSL::gsl GSL::gslcblas)

pybind11_add_module(_khuri_omnes
    "${SOURCE_DIR}/cauchy.cpp"
    "${SOURCE_DIR}/gsl_interface.cpp"
    "${BINDING_DIR}/omnes_bindings.cpp")
target_link_libraries(_khuri_omnes PRIVATE GSL::gsl GSL::gslcblas)

pybind11_add_module(_khuri_khuri_treiman
    "${SOURCE_DIR}/cauchy.cpp"
    "${SOURCE_DIR}/curved_omnes.cpp"
    "${SOURCE_DIR}/grid.cpp"
    "${SOURCE_DIR}/gsl_interface.cpp"
    "${SOURCE_DIR}/kernel.cpp"
    "${SOURCE_DIR}/partial_wave.cpp"
    "${SOURCE_DIR}/piecewise.cpp"
    "${BINDING_DIR}/khuri_treiman_bindings.cpp")
target_link_libraries(_khuri_khuri_treiman PRIVATE GSL::gsl GSL::gslcblas)

pybind11_add_module(_khuri_pinocchio
    "${SOURCE_DIR}/gsl_interface.cpp"
    "${SOURCE_DIR}/cauchy.cpp"
    "${SOURCE_DIR}/phase.cpp"
    "${SOURCE_DIR}/asymptotic.cpp"
    "${SOURCE_DIR}/path.cpp"
    "${SOURCE_DIR}/array.cpp"
    "${SOURCE_DIR}/matching.cpp"
    "${SOURCE_DIR}/splined_omnes.cpp"
    "${SOURCE_DIR}/splined_path.cpp"
    "${SOURCE_DIR}/path_eta_pi_pi.cpp"
    "${SOURCE_DIR}/angular_average.cpp"
    "${SOURCE_DIR}/dispersion_integral.cpp"
    "${SOURCE_DIR}/basis_function.cpp"
    "${SOURCE_DIR}/iterative_solution.cpp"
    "${BINDING_DIR}/khuri_pinocchio_bindings.cpp")
target_link_libraries(_khuri_pinocchio PRIVATE GSL::gsl GSL::gslcblas)

pybind11_add_module(_khuri_chpt
    "${SOURCE_DIR}/chpt.cpp"
    "${BINDING_DIR}/chpt_bindings.cpp")

pybind11_add_module(_khuri_iam
    "${SOURCE_DIR}/chpt.cpp"
    "${SOURCE_DIR}/iam.cpp"
    "${BINDING_DIR}/iam_bindings.cpp")

pybind11_add_module(_khuri_mandelstam
    "${BINDING_DIR}/mandelstam_bindings.cpp")

pybind11_add_module(_khuri_curved_omnes
    "${SOURCE_DIR}/cauchy.cpp"
    "${SOURCE_DIR}/curved_omnes.cpp"
    "${SOURCE_DIR}/grid.cpp"
    "${SOURCE_DIR}/gsl_interface.cpp"
    "${SOURCE_DIR}/piecewise.cpp"
    "${BINDING_DIR}/curved_omnes_bindings.cpp")
target_link_libraries(_khuri_curved_omnes PRIVATE GSL::gsl GSL::gslcblas)

pybind11_add_module(_khuri_path
    "${SOURCE_DIR}/cauchy.cpp"
    "${SOURCE_DIR}/gsl_interface.cpp"
    "${SOURCE_DIR}/path.cpp"
    "${BINDING_DIR}/path_bindings.cpp")
target_link_libraries(_khuri_path PRIVATE GSL::gsl GSL::gslcblas)
