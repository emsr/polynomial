/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../include -o test_jacobi_roots test_jacobi_roots.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_jacobi_roots > test_jacobi_roots.txt
*/

#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <ext/solver_jenkins_traub.h>
#include <ext/polynomial.h>

// Just grab this Jacobi polynomial from tr29124_test.
namespace std
{
namespace __detail
{
  /**
   * Return the Jacobi polynomial as a polynomial.
   *
   * @tparam _Tp The real type of the argument and degree parameters.
   * @param[in]  n  The order of the Jacobi polynomial
   * @param[in]  alpha1  The first parameter of the Jacobi polynomial
   * @param[in]  beta1  The second parameter of the Jacobi polynomial
   */
  template<typename _Tp>
    __gnu_cxx::_Polynomial<_Tp>
    __jacobi_poly(unsigned int __n, _Tp __alpha1, _Tp __beta1)
    {
      __gnu_cxx::_Polynomial<_Tp> __poly;
      const __gnu_cxx::_Polynomial<_Tp> __arg{_Tp{0.5L}, _Tp{-0.5L}};

      if (std::isnan(__alpha1) || std::isnan(__beta1))
	return __poly;

      auto __term = __gnu_cxx::_Polynomial<_Tp>{1};
      __poly += __term;
      if (__n == 0)
	return __poly;

      auto __fact = _Tp{1};
      const auto __ab = __alpha1 + __beta1;
      for (unsigned int __k = 1; __k <= __n; ++__k)
	{
	  __fact *= _Tp(__alpha1 + __k) / _Tp(__k);
	  __term *= _Tp(-__n + __k - 1) / _Tp(__k)
		  * _Tp(__n + __k + __ab) / _Tp(__alpha1 + __k)
		  * __arg;
	  __poly += __term;
	}

      return __fact * __poly;
    }
} // namespace std
} // namespace __detail

template<typename _Tp>
  void
  test_jacobi_roots(unsigned n, _Tp alpha1, _Tp beta1, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<_Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; alpha = " << alpha1
	      << "; beta = " << beta1 << '\n';

    const auto poly = std::__detail::__jacobi_poly(n, alpha1, beta1);
    auto coef = poly.coefficients();
    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';

    std::reverse(coef.begin(), coef.end());

    auto jt = __gnu_cxx::_JenkinsTraubSolver(coef);
    auto roots = jt.solve();
    std::cout << "\nThe roots are:\n";
    for (const auto& z : roots)
      {
	if (z.index() == 0)
	  continue;
	else if (z.index() == 1)
	  std::cout << ' ' << std::setw(w) << std::get<1>(z)
		    << ' ' << std::setw(w) << 0.0
		    << ' ' << std::setw(w) << poly(std::get<1>(z))
		    << '\n';
	else
	  std::cout << ' ' << std::setw(w) << std::real(std::get<2>(z))
		    << ' ' << std::setw(w) << std::imag(std::get<2>(z))
		    << ' ' << std::setw(w) << poly(std::get<2>(z))
		    << '\n';
      }

    gp << std::setprecision(prec);
    gp << "\n\n";
    gp << "# n = " << n
       << "; alpha = " << alpha1
       << "; beta = " << beta1 << '\n';
    for (const auto& z : roots)
      {
	if (z.index() == 0)
	  continue;
	else if (z.index() == 1)
	  gp << ' ' << std::setw(w) << std::get<1>(z)
	     << ' ' << std::setw(w) << 0.0
	     << '\n';
	else
	  gp << ' ' << std::setw(w) << std::real(std::get<2>(z))
	     << ' ' << std::setw(w) << std::imag(std::get<2>(z))
	     << '\n';
      }
    gp << std::flush;
  }

/**
 * Numerical Methods for Special Functions, Gil, Segura, Temme, pp. 192.
 */
template<typename _Tp>
  void
  run()
  {
    std::ofstream gp("roots.gp");

    unsigned n = 50;
    _Tp alpha1, beta1;

    alpha1 = _Tp{2.0L};
    beta1 = _Tp{-42.5L};
    test_jacobi_roots(n, alpha1, beta1, gp);

    alpha1 = _Tp{2.0L};
    beta1 = _Tp{-52.1L};
    test_jacobi_roots(n, alpha1, beta1, gp);

    alpha1 = _Tp{2.0L};
    beta1 = _Tp{-63.5L};
    test_jacobi_roots(n, alpha1, beta1, gp);

    // Flip alpha and beta.

    alpha1 = _Tp{2.0L};
    beta1 = _Tp{-42.5L};
    test_jacobi_roots(n, beta1, alpha1, gp);

    alpha1 = _Tp{2.0L};
    beta1 = _Tp{-52.1L};
    test_jacobi_roots(n, beta1, alpha1, gp);

    alpha1 = _Tp{2.0L};
    beta1 = _Tp{-63.5L};
    test_jacobi_roots(n, beta1, alpha1, gp);
  }

int
main()
{
  run<long double>();
}
