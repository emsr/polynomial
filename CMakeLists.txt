cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_polynomial
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

add_custom_target(make_cxx_polynomial_output_dir ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory output)

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs_cxx_polynomial
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_polynomial" VERBATIM
  )
endif (DOXYGEN_FOUND)

# Header-only library.
add_library(cxx_polynomial INTERFACE)
target_include_directories(cxx_polynomial INTERFACE include)

add_executable(test_bairstow test_bairstow.cpp)
target_link_libraries(test_bairstow cxx_polynomial quadmath)

add_executable(test_horner test_horner.cpp)
target_link_libraries(test_horner cxx_polynomial quadmath)

add_executable(test_jenkins_traub test_jenkins_traub.cpp)
target_link_libraries(test_jenkins_traub cxx_polynomial quadmath)

add_executable(test_polynomial test_polynomial.cpp)
target_link_libraries(test_polynomial cxx_polynomial quadmath)

add_executable(test_polynomial_root test_polynomial_root.cpp)
target_link_libraries(test_polynomial_root cxx_polynomial quadmath)

add_executable(test_rational_polynomial test_rational_polynomial.cpp)
target_link_libraries(test_rational_polynomial cxx_polynomial quadmath)

add_executable(test_solution test_solution.cpp)
target_link_libraries(test_solution cxx_polynomial quadmath)

add_executable(test_solvers test_solvers.cpp)
target_link_libraries(test_solvers cxx_polynomial quadmath)

add_executable(test_static_polynomial test_static_polynomial.cpp)
target_link_libraries(test_static_polynomial cxx_polynomial quadmath)

add_executable(test_laguerre_step test_laguerre_step.cpp)
target_link_libraries(test_laguerre_step cxx_polynomial quadmath)

add_executable(test_quadratic_step test_quadratic_step.cpp)
target_link_libraries(test_quadratic_step cxx_polynomial quadmath)

# Requires tr29124...

if (FOUND_TR29124)
  add_executable(test_jacobi_roots test_jacobi_roots.cpp)
  target_include_directories(test_jacobi_roots PRIVATE
    ../include
    ../cxx_math_constants/include
    ../cxx_fp_utils/include)
  target_link_libraries(test_jacobi_roots cxx_polynomial quadmath)
endif(FOUND_TR29124)
