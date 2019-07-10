cmake_minimum_required (VERSION 3.10)

include(CheckCXXCompilerFlag)

project(
  cxx_polynomial
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS -Wno-psabi)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target(docs_polynomial
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation for cxx_polynomial" VERBATIM
  )
endif (DOXYGEN_FOUND)

add_executable(test_bairstow test_bairstow.cpp)
target_include_directories(test_bairstow PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_bairstow quadmath)

add_executable(test_horner test_horner.cpp)
target_include_directories(test_horner PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_horner quadmath)

add_executable(test_jenkins_traub test_jenkins_traub.cpp)
target_include_directories(test_jenkins_traub PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_jenkins_traub quadmath)

add_executable(test_polynomial test_polynomial.cpp)
target_include_directories(test_polynomial PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_polynomial quadmath)

add_executable(test_polynomial_root test_polynomial_root.cpp)
target_include_directories(test_polynomial_root PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_polynomial_root quadmath)

add_executable(test_rational_polynomial test_rational_polynomial.cpp)
target_include_directories(test_rational_polynomial PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_rational_polynomial quadmath)

add_executable(test_solution test_solution.cpp)
target_include_directories(test_solution PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_solution quadmath)

add_executable(test_solvers test_solvers.cpp)
target_include_directories(test_solvers PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_solvers quadmath)

add_executable(test_static_polynomial test_static_polynomial.cpp)
target_include_directories(test_static_polynomial PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_static_polynomial quadmath)

add_executable(test_laguerre_step test_laguerre_step.cpp)
target_include_directories(test_laguerre_step PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_laguerre_step quadmath)

add_executable(test_quadratic_step test_quadratic_step.cpp)
target_include_directories(test_quadratic_step PRIVATE include ../cxx_math_constants/include)
target_link_libraries(test_quadratic_step quadmath)

# Requires tr29124...

if (FOUND_TR29124)
  add_executable(test_jacobi_roots test_jacobi_roots.cpp)
  target_include_directories(test_jacobi_roots PRIVATE
    include
    ../include
    ../cxx_math_constants/include
    ../cxx_fp_utils/include)
  target_link_libraries(test_jacobi_roots quadmath)
endif(FOUND_TR29124)
