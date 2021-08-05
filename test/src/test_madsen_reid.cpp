
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

    a1.resize(m + 2);

    for (int k = 1; k <= m + 1; ++k)
      {
        std::cout << "Enter coefficient " << k - 1 << ": ";
        std::cin >> a1[m + 2 - k];
      }

    SolverMadsenReid<Real> smr(a1);
    auto root = smr.solve();

    std::cout << '\n';
    for (int k = 1; k <= m; ++k)
      {
        std::cout << "Zero " << k - 1 << ": " << root[k] << '\n';
      }

    int num_errors = 0;
    const auto tol = std::sqrt(std::numeric_limits<Real>::epsilon());

    std::cout << "\nSolution tests:\n";
    // Remember to reverse the polynomial coefficients!
    __gnu_cxx::_Polynomial<Cmplx> poly(a1.rbegin(), a1.rend() - 1);
    std::cout << "P(z) = " << poly << '\n';
    //for (const auto& z : zeros)
    for (int k = 1; k <= m; ++k)
      {
        const auto z = root[k];
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
