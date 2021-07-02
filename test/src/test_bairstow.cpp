
#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>

#include "ext/polynomial.h"
#include "ext/solver_bairstow.h"

template<typename _Real>
  int
  test_bairstow()
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);
    //const auto w = 6 + std::cout.precision();
    //const auto cw = 4 + 2 * w;

    int order = 0;
    int MAX_TERMS = 1000;
    while ((order < 2) || (order > MAX_TERMS - 1))
      {
	std::cout << "Polynomial order (2 - " << MAX_TERMS - 1 << "; 0 to quit): ";
	std::cin >> order;
        if (order <= 0)
          break;
      }
    std::vector<_Real> a(order + 1);

    std::cout << "Enter coefficients, high order to low order.\n";
    for (int i = 0; i <= order; ++i)
      {
	std::cout << "a[" << (order - i) << "] = ";
	std::cin >> a[i];
      }
    std::cout << "\nP : ( ";
    for (int i = 0; i < order; ++i)
      std::cout << a[order - i] << ", ";
    std::cout << a[0] << " )\n";

    __gnu_cxx::_BairstowSolver bairstow(a, 123456);
    const auto zeros = bairstow.solve();
    std::cout << "\nThe zeros are:\n";
    for (const auto& z : zeros)
      std::cout << z << '\n';

    const auto eq = bairstow.equations();
    std::cout << "\nThe quadratic factors are:\n";
    for (unsigned int p = 0; p < eq.size() / 2; ++p)
      std::cout << "t^2 + " << eq[2 * p + 1] << " t + " << eq[2 * p] << '\n';
    if ((eq.size() % 2) == 1)
      std::cout << "The linear term is: \nt - " << eq.back() << '\n';

    int num_errors = 0;
    const auto tol = std::sqrt(std::numeric_limits<_Real>::epsilon());

    std::cout << "\nSolution tests:\n";
    __gnu_cxx::_Polynomial<_Real> poly(a.begin(), a.end());
    for (const auto& z : zeros)
      {
	const auto idx = z.index();
	std::cout << "P(" << z << ") = ";
	if (idx == 1)
	  {
	    const auto P = poly(std::get<1>(z));
	    std::cout << P << '\n';
	    if (const auto aP = std::abs(P); aP > tol)
	      {
		++num_errors;
		std::cout << "Fail: |P| = " << aP << '\n';
	      }
          }
	else if (idx == 2)
          {
	    const auto P = poly(std::get<2>(z));
	    std::cout << P << '\n';
	    if (const auto aP = std::abs(P); aP > tol)
	      {
		++num_errors;
		std::cout << "Fail: |P| = " << aP << '\n';
	      }
          }
      }

    return num_errors;
  }

int
main()
{
  int num_errors = 0;

  std::cout << "\ndouble\n======\n";
  num_errors += test_bairstow<double>();

  std::cout << "\nlong double\n===========\n";
  num_errors += test_bairstow<long double>();

  std::cout << "\nfloat\n=====\n";
  num_errors += test_bairstow<float>();

  return num_errors;
}
