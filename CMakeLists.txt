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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/complex_solver1.in ${CMAKE_CURRENT_BINARY_DIR}/complex_solver1.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/complex_solver2.in ${CMAKE_CURRENT_BINARY_DIR}/complex_solver2.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver1.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver1.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver2.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver2.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver3.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver3.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver4.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver4.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver5.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver5.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver6.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver6.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver7.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver7.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver8.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver8.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_solver9.in ${CMAKE_CURRENT_BINARY_DIR}/test_solver9.in COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/input/test_bairstow.in ${CMAKE_CURRENT_BINARY_DIR}/test_bairstow.in COPYONLY)

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

add_executable(test_bairstow test/src/test_bairstow.cpp)
target_link_libraries(test_bairstow cxx_polynomial quadmath)
add_test(NAME run_test_bairstow COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_bairstow < ${CMAKE_CURRENT_BINARY_DIR}/test_bairstow.in > output/test_bairstow.txt")

add_executable(test_horner test/src/test_horner.cpp)
target_link_libraries(test_horner cxx_polynomial quadmath)
add_test(NAME run_test_horner COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_horner > output/test_horner.txt")

add_executable(test_jenkins_traub test/src/test_jenkins_traub.cpp)
target_link_libraries(test_jenkins_traub cxx_polynomial quadmath)
add_test(NAME run_test_jenkins_traub_1 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver1.in > output/test_jenkins_traub_1.txt")
add_test(NAME run_test_jenkins_traub_2 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver2.in > output/test_jenkins_traub_2.txt")
add_test(NAME run_test_jenkins_traub_3 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver3.in > output/test_jenkins_traub_3.txt")
add_test(NAME run_test_jenkins_traub_4 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver4.in > output/test_jenkins_traub_4.txt")
add_test(NAME run_test_jenkins_traub_5 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver5.in > output/test_jenkins_traub_5.txt")
add_test(NAME run_test_jenkins_traub_6 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver6.in > output/test_jenkins_traub_6.txt")
add_test(NAME run_test_jenkins_traub_7 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver7.in > output/test_jenkins_traub_7.txt")
add_test(NAME run_test_jenkins_traub_8 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver8.in > output/test_jenkins_traub_8.txt")
add_test(NAME run_test_jenkins_traub_9 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub < ${CMAKE_CURRENT_BINARY_DIR}/test_solver9.in > output/test_jenkins_traub_9.txt")

add_executable(test_jenkins_traub_complex test/src/test_jenkins_traub_complex.cpp)
target_link_libraries(test_jenkins_traub_complex cxx_polynomial quadmath)
add_test(NAME run_test_jenkins_traub_complex_1 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver1.in > output/test_jenkins_traub_complex_1.txt")
add_test(NAME run_test_jenkins_traub_complex_2 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver2.in > output/test_jenkins_traub_complex_2.txt")
add_test(NAME run_test_jenkins_traub_complex_3 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver3.in > output/test_jenkins_traub_complex_3.txt")
add_test(NAME run_test_jenkins_traub_complex_4 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver4.in > output/test_jenkins_traub_complex_4.txt")
add_test(NAME run_test_jenkins_traub_complex_5 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver5.in > output/test_jenkins_traub_complex_5.txt")
add_test(NAME run_test_jenkins_traub_complex_6 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver6.in > output/test_jenkins_traub_complex_6.txt")
add_test(NAME run_test_jenkins_traub_complex_7 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver7.in > output/test_jenkins_traub_complex_7.txt")
add_test(NAME run_test_jenkins_traub_complex_8 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver8.in > output/test_jenkins_traub_complex_8.txt")
add_test(NAME run_test_jenkins_traub_complex_9 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/test_solver9.in > output/test_jenkins_traub_complex_9.txt")
add_test(NAME test_jenkins_traub_complex_complex_1 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/complex_solver1.in > output/test_jenkins_traub_complex_complex_1.txt")
add_test(NAME test_jenkins_traub_complex_complex_2 COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_jenkins_traub_complex < ${CMAKE_CURRENT_BINARY_DIR}/complex_solver2.in > output/test_jenkins_traub_complex_complex_2.txt")

add_executable(test_polynomial test/src/test_polynomial.cpp)
target_link_libraries(test_polynomial cxx_polynomial quadmath)
add_test(NAME run_test_polynomial COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_polynomial > output/test_polynomial.txt")

add_executable(test_polynomial_root test/src/test_polynomial_root.cpp)
target_link_libraries(test_polynomial_root cxx_polynomial quadmath)
add_test(NAME run_test_polynomial_root COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_polynomial_root > output/test_polynomial_root.txt")

add_executable(test_rational_polynomial test/src/test_rational_polynomial.cpp)
target_link_libraries(test_rational_polynomial cxx_polynomial quadmath)
add_test(NAME run_test_rational_polynomial COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_rational_polynomial > output/test_rational_polynomial.txt")

add_executable(test_solution test/src/test_solution.cpp)
target_link_libraries(test_solution cxx_polynomial quadmath)

add_executable(test_solvers test/src/test_solvers.cpp)
target_link_libraries(test_solvers cxx_polynomial quadmath)
add_test(NAME run_test_solvers COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_solvers > output/test_solvers.txt")

add_executable(test_static_polynomial test/src/test_static_polynomial.cpp)
target_link_libraries(test_static_polynomial cxx_polynomial quadmath)
add_test(NAME run_test_static_polynomial COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_static_polynomial > output/test_static_polynomial.txt")

add_executable(test_laguerre_step test/src/test_laguerre_step.cpp)
target_link_libraries(test_laguerre_step cxx_polynomial quadmath)
add_test(NAME run_test_laguerre_step COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_laguerre_step > output/test_laguerre_step.txt")

add_executable(test_quadratic_step test/src/test_quadratic_step.cpp)
target_link_libraries(test_quadratic_step cxx_polynomial quadmath)
add_test(NAME run_test_quadratic_step COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_quadratic_step > output/test_quadratic_step.txt")

add_executable(test_derivative_integral test/src/test_derivative_integral.cpp)
target_link_libraries(test_derivative_integral cxx_polynomial quadmath)
add_test(NAME run_test_derivative_integral COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_derivative_integral > output/test_derivative_integral.txt")

add_executable(test_shift test/src/test_shift.cpp)
target_link_libraries(test_shift cxx_polynomial quadmath)
add_test(NAME run_test_shift COMMAND bash -c "${CMAKE_BINARY_DIR}/bin/test_shift > output/test_shift.txt")

# Requires tr29124...

if (FOUND_TR29124)
  add_executable(test_jacobi_roots test_jacobi_roots.cpp)
  target_include_directories(test_jacobi_roots PRIVATE
    ../include
    ../cxx_math_constants/include
    ../cxx_fp_utils/include)
  target_link_libraries(test_jacobi_roots cxx_polynomial quadmath)
endif(FOUND_TR29124)
