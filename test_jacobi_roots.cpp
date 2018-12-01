/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -I../include -o test_jacobi_roots test_jacobi_roots.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_jacobi_roots > test_jacobi_roots.txt
*/

#include <vector>
#include <algorithm>
#include <fstream>

#include <ext/solver_jenkins_traub.h>
#include <ext/math_const.h>
#include <bits/specfun.h>

template<typename _Tp>
  void
  test(unsigned n, _Tp alpha1, _Tp beta1, std::ofstream& gp)
  {
    auto poly = std::__detail::__jacobi_poly(n, alpha1, beta1);
    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : poly)
      std::cout << c << '\n';

    std::reverse(poly.begin(), poly.end());

    auto jt = __gnu_cxx::_JenkinsTraubSolver(poly);
    auto roots = jt.solve();
    std::cout << "\nThe roots are:\n";
    for (const auto& z : roots)
      std::cout << z << '\n';

    gp << "\n\n";
    for (const auto& z : roots)
      {
	if (z.index() == 0)
	  continue;
	else if (z.index() == 1)
	  gp << "  " << std::get<1>(z)
	     << "  " << 0.0 << '\n';
	else
	  gp << "  " << std::real(std::get<2>(z))
	     << "  " << std::imag(std::get<2>(z)) << '\n';
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

    alpha1 = 2.0L;
    beta1 = -42.5L;
    test(n, alpha1, beta1, gp);

    alpha1 = 2.0L;
    beta1 = -52.1L;
    test(n, alpha1, beta1, gp);

    alpha1 = 2.0L;
    beta1 = -63.5L;
    test(n, alpha1, beta1, gp);

    alpha1 = 2.0L;
    beta1 = -42.5L;
    test(n, beta1, alpha1, gp);

    alpha1 = 2.0L;
    beta1 = -52.1L;
    test(n, beta1, alpha1, gp);

    alpha1 = 2.0L;
    beta1 = -63.5L;
    test(n, beta1, alpha1, gp);
  }

int
main()
{
  run<long double>();
}
