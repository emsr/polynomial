/*
g++ -std=c++20 -I. -Iinclude -o test_solver_skowron_gould test_solver_skowron_gould.cpp
*/

#include <iostream>
#include <complex>

#include <emsr/polynomial.h>
#include <emsr/notsospecfun.h>
#include <emsr/solver_low_degree.h>
#include <emsr/solver_skowron_gould.h>

template<typename Real>
  int
  test_skowron_gould()
  {
    using Cmplx = std::complex<Real>;

    std::cout.precision(std::numeric_limits<Real>::digits10);

    int m;
    std::cout << "Enter degree: ";
    std::cin >> m;
    if (m <= 0)
      return 0;

    std::vector<Cmplx> coeff(m + 1);

    for (int k = 0; k <= m; ++k)
      {
        std::cout << "Enter coefficient " << k << ": ";
        std::cin >> coeff[k];
      }

    //SolverSkowronGould<Real> ssg(coeff);
    //auto root = ssg.solve();
    std::vector<Cmplx> root(m);
    bool polish_roots_after = true;
    bool use_roots_as_starting_points = false;
    cmplx_roots_gen(root, coeff, m, polish_roots_after, use_roots_as_starting_points);

    std::cout << '\n';
    for (unsigned int k = 0; k < root.size(); ++k)
      std::cout << "Zero " << k << ": " << root[k] << '\n';

    int num_errors = 0;
    const auto tol = std::sqrt(std::numeric_limits<Real>::epsilon());

    std::cout << "\nSolution tests:\n";
    // Remember to reverse the polynomial coefficients!
    emsr::Polynomial<Cmplx> poly(coeff.begin(), coeff.end());
    std::cout << "P(z) = " << poly << '\n';
    for (const auto& z : root)
      {
	std::cout << "P(" << z << ") = ";
	const auto P = poly(z);
	std::cout << P << '\n';
	if (const auto aP = std::abs(P); aP > tol)
	  {
	    ++num_errors;
	    std::cout << "Fail: |P| = " << aP << '\n';
	  }
      }

    return num_errors;
  }

int
main()
{
  int num_errors = 0;

  std::cout << "\ndouble\n======\n" << std::flush;
  num_errors += test_skowron_gould<double>();

  std::cout << "\nlong double\n===========\n" << std::flush;
  num_errors += test_skowron_gould<long double>();

  std::cout << "\nfloat\n=====\n" << std::flush;
  num_errors += test_skowron_gould<float>();

  return num_errors;
}
