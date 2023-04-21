cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# BEGIN GOOGLETEST SETUP ######################################################

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# END GOOGLETEST SETUP ########################################################

add_definitions(-Wall -pedantic)

set(SOURCE_DIR "../src")
set(INCLUDE_DIR "../include")

include_directories(${INCLUDE_DIR})

add_executable(test_gsl
    test_gsl_interface.cpp
    "${SOURCE_DIR}/gsl_interface.cpp"
    "${SOURCE_DIR}/iam.cpp"
    "${SOURCE_DIR}/asymptotic.cpp"
    "${SOURCE_DIR}/chpt.cpp")
target_link_libraries(test_gsl PRIVATE gtest Threads::Threads GSL::gsl GSL::gslcblas)


add_executable(test_cauchy
    test_cauchy.cpp
    ${SOURCE_DIR}/gsl_interface.cpp
    ${SOURCE_DIR}/cauchy.cpp)
target_link_libraries(test_cauchy PRIVATE gmock gtest Threads::Threads GSL::gsl GSL::gslcblas)

enable_testing()

add_test(NAME gsl
    COMMAND ./test_gsl)

add_test(NAME cauchy
    COMMAND ./test_cauchy)
