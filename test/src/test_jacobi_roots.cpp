
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <cmath> // FIXME: For isnan for math_util.h
#include <emsr/math_util.h>
#include <emsr/solver_jenkins_traub.h>
#include <emsr/polynomial.h>

#include <bits/float128_io.h>
#include <bits/float128_math.h>

namespace std
{
namespace __detail
{
  /**
   * Return the Jacobi polynomial as a polynomial.
   *
   * @tparam Tp The real type of the argument and degree parameters.
   * @param[in]  n  The degree of the Jacobi polynomial
   * @param[in]  alpha1  The first parameter of the Jacobi polynomial
   * @param[in]  beta1  The second parameter of the Jacobi polynomial
   */
  template<typename Tp>
    emsr::Polynomial<Tp>
    jacobi_poly(unsigned int n, Tp alpha1, Tp beta1)
    {
      emsr::Polynomial<Tp> poly;
      const emsr::Polynomial<Tp> arg({Tp{0.5L}, Tp{-0.5L}});

      if (std::isnan(alpha1) || std::isnan(beta1))
	return poly;

      auto term = emsr::Polynomial<Tp>{1};
      poly += term;
      if (n == 0)
	return poly;

      auto fact = Tp{1};
      const auto ab = alpha1 + beta1;

      int m = int(n);
      const auto Maybe = emsr::fp_is_integer(n + ab);
      if (Maybe && Maybe() < 0 && -Maybe() < m)
	m = -Maybe();

      for (int k = 1; k <= m; ++k)
	{
	  fact *= Tp(alpha1 + k) / Tp(k);

	  term *= (Tp(-m + k - 1) / Tp(k))
		  * (Tp(m + k + ab) / Tp(alpha1 + k))
		  * arg;

	  poly += term;
	}

      return fact * poly;
    }

  /**
   * Highest degree term coefficient.
   */
  template<typename Tp>
    Tp
    jacobi_norm(unsigned int n, Tp alpha1, Tp beta1)
    {
      int sgam1, sgam2;
      const auto lgam1 = lgamma_r(Tp(2 * n + alpha1 + beta1 + 1), &sgam1);
      const auto lgam2 = lgamma_r(Tp(n + alpha1 + beta1 + 1), &sgam2);
      return sgam1 * sgam2 * std::exp(lgam1 - std::lgamma(Tp(n + 1))
   				    - lgam2 - Tp(n) * std::log(Tp{2}));
    }

} // namespace __detail
} // namespace std

template<typename Tp>
  void
  test_jacobi_roots(unsigned n, Tp alpha1, Tp beta1, std::ofstream& gp)
  {
    const auto prec = std::numeric_limits<Tp>::digits10;
    const auto w = 6 + prec;

    std::cout << std::setprecision(prec);
    std::cout << "\n\n";
    std::cout << " n = " << n
	      << "; alpha = " << alpha1
	      << "; beta = " << beta1 << '\n';

    const auto poly = std::detail::jacobi_poly(n, alpha1, beta1);
    auto coef = poly.coefficients();
    std::cout << "\nThe polynomial coefficients are:\n";
    for (const auto& c : coef)
      std::cout << std::setw(w) << c << '\n';
    std::cout << "\nMax coefficient: " << std::__detail::jacobi_norm(n, alpha1, beta1) << '\n';

    std::reverse(coef.begin(), coef.end());

    auto jt = emsr::JenkinsTraubSolver(coef);
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
template<typename Tp>
  void
  run()
  {
    std::ofstream gp("jacobi_roots.dat");

    unsigned n = 50;
    Tp alpha1, beta1;

    alpha1 = Tp{2.0L};
    beta1 = Tp{-42.5L};
    test_jacobi_roots(n, alpha1, beta1, gp);

    alpha1 = Tp{2.0L};
    beta1 = Tp{-52.0L};
    test_jacobi_roots(n, alpha1, beta1, gp);

    alpha1 = Tp{2.0L};
    beta1 = Tp{-63.5L};
    test_jacobi_roots(n, alpha1, beta1, gp);

    // Flip alpha and beta.

    alpha1 = Tp{2.0L};
    beta1 = Tp{-42.5L};
    test_jacobi_roots(n, beta1, alpha1, gp);

    alpha1 = Tp{2.0L};
    beta1 = Tp{-52.0L};
    test_jacobi_roots(n, beta1, alpha1, gp);

    alpha1 = Tp{2.0L};
    beta1 = Tp{-63.5L};
    test_jacobi_roots(n, beta1, alpha1, gp);
  }

int
main()
{
  //run<long double>();

  run<__float128>();
}
