
#include <vector>
#include <limits>
#include <complex>
#include <iostream>
#include <iomanip>

#include <ext/polynomial.h>
#include <ext/solver_madsen_reid.h>

template<typename Real>
  int
  test_madsen_reid()
  {
    using Cmplx = std::complex<Real>;

    std::cout.precision(std::numeric_limits<Real>::digits10);

    int m;
    std::vector<Cmplx> a1;

    std::cout << "Enter degree: ";
    std::cin >> m;
    if (m <= 0)
      return 0;

    a1.resize(m + 1);

    for (int k = 0; k <= m; ++k)
      {
        std::cout << "Enter coefficient " << k << ": ";
        std::cin >> a1[m - k];
      }

    SolverMadsenReid<Real> smr(a1);
    auto root = smr.solve();

    std::cout << '\n';
    for (unsigned int k = 0; k < root.size(); ++k)
      std::cout << "Zero " << k << ": " << root[k] << '\n';

    int num_errors = 0;
    const auto tol = std::sqrt(std::numeric_limits<Real>::epsilon());

    std::cout << "\nSolution tests:\n";
    // Remember to reverse the polynomial coefficients!
    emsr::Polynomial<Cmplx> poly(a1.begin(), a1.end());
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
  num_errors += test_madsen_reid<double>();

  std::cout << "\nlong double\n===========\n" << std::flush;
  num_errors += test_madsen_reid<long double>();

  std::cout << "\nfloat\n=====\n" << std::flush;
  num_errors += test_madsen_reid<float>();

  return num_errors;
}
