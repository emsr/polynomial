/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bairstow test_bairstow.cpp
./test_bairstow < test_bairstow.in > test_bairstow.txt
*/

#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>

//#include "bits/specfun_state.h"
#include "ext/polynomial.h"
#include "ext/solver_low_degree.h"

namespace __gnu_cxx
{

///
/// @todo: If you *don't* reverse the input array you solve for 1/z.
/// @todo Take a max_error.  If this->_M_eps grows larger than max_error throw.
///
template<typename _Real>
  class _BairstowSolver
  {
  public:

    _BairstowSolver(const std::vector<_Real>& __coeff,
		    unsigned int __seed = std::random_device()())
    : _M_coeff{__coeff.rbegin(), __coeff.rend()},
      _M_b(__coeff.size()), _M_c(__coeff.size()),
      _M_order(__coeff.size() - 1),
      _M_urng(__seed), _M_pdf(_Real{0}, _Real{2})
    {
      if (this->_M_coeff.size() == 0)
	std::__throw_domain_error("_BairstowSolver: "
				  "Coefficient size must be nonzero.");

      if (this->_M_coeff[0] == _Real{0})
	std::__throw_domain_error("_BairstowSolver: "
				  "Leading-order coefficient must be nonzero.");

      const auto __scale = this->_M_coeff[0];
      for (int __i = 0; __i <= this->_M_order; ++__i)
	this->_M_coeff[__i] /= __scale;

      this->_M_zero.reserve(this->_M_coeff.size());
    }

    std::vector<solution_t<_Real>> solve();
    std::vector<_Real> equations() const;

  private:

    void _M_iterate();

    template<int _Index>
    static void _S_refine_quadratic_eqn(_Real& __dr, _Real& __r,
					_Real& __ds, _Real& __s,
					std::array<solution_t<_Real>, 2>& __w);

    void
    _M_add_zero(solution_t<_Real> __z)
    {
      this->_M_zero.push_back(__z);
      --this->_M_order;
    }

    static constexpr _Real _S_ratio
	= _Real{std::numeric_limits<_Real>::digits}
	/ std::numeric_limits<double>::digits;
    static constexpr int _S_max_rand_iter = 200 * _S_ratio;
    static constexpr int _S_max_error_iter = 500 * _S_ratio;
    static constexpr auto _S_eps_factor = _Real{10} * _S_ratio;
    static constexpr auto _S_eps
      = _S_eps_factor * std::numeric_limits<_Real>::epsilon();

    std::vector<_Real> _M_coeff;
    std::vector<_Real> _M_b;
    std::vector<_Real> _M_c;
    std::vector<solution_t<_Real>> _M_zero;
    _Real _M_eps = _S_eps;
    int _M_order;
    bool _M_precision_error = false;
    std::mt19937 _M_urng;
    std::uniform_real_distribution<_Real> _M_pdf;
  };

/**
 * Attempt to find a quadratic factor of the polynomial by Bairstow's method.
 */
template<typename _Real>
  void
  _BairstowSolver<_Real>::_M_iterate()
  {
    auto __r = _Real{0};
    auto __s = _Real{0};
    auto __dr = _Real{1};
    auto __ds = _Real{0};

    int __iter = 1;
    this->_M_precision_error = false;

    while (std::abs(__dr) + std::abs(__ds) > this->_M_eps)
      {
	if (__iter % _S_max_rand_iter == 0)
	  __r = this->_M_pdf(this->_M_urng);
	if (__iter % _S_max_error_iter == 0)
	  {
	    this->_M_eps *= _S_eps_factor;
	    this->_M_precision_error = true;
	    std::cout << "Loss of precision: " << this->_M_eps << '\n';
	  }

	this->_M_b[1] = this->_M_coeff[1] - __r;
	this->_M_c[1] = this->_M_b[1] - __r;
	for (int __i = 2; __i <= this->_M_order; ++__i)
	  {
	    this->_M_b[__i] = this->_M_coeff[__i]
			- __r * this->_M_b[__i - 1] - __s * this->_M_b[__i - 2];
	    this->_M_c[__i] = this->_M_b[__i]
			- __r * this->_M_c[__i - 1] - __s * this->_M_c[__i - 2];
	  }

	auto __dn = this->_M_c[this->_M_order - 1]
		  * this->_M_c[this->_M_order - 3]
		  - this->_M_c[this->_M_order - 2]
		  * this->_M_c[this->_M_order - 2];
	if (std::abs(__dn) < std::numeric_limits<_Real>::epsilon())
	  {
	    __dr = _Real{1};
	    __ds = _Real{1};
	  }
	else
	  {
	    auto __drn = this->_M_b[this->_M_order]
		       * this->_M_c[this->_M_order - 3]
		       - this->_M_b[this->_M_order - 1]
		       * this->_M_c[this->_M_order - 2];
	    auto __dsn = this->_M_b[this->_M_order - 1]
		       * this->_M_c[this->_M_order - 1]
		       - this->_M_b[this->_M_order]
		       * this->_M_c[this->_M_order - 2];
	    __dr = __drn / __dn;
	    __ds = __dsn / __dn;
	  }

	__r += __dr;
	__s += __ds;
	++__iter;
	if (std::abs(__dr) + std::abs(__ds) <= this->_M_eps)
	  {
	    // Before exiting, solve the quadratic, refine the roots,
	    // multiply out the resulting polynomial, extract
	    // the new coefficients, get new dr, r, ds, s.
	    auto __w = __quadratic(__s, __r, _Real{1});
	    if (__w[0].index() == 1 && __w[1].index() == 1)
	      _S_refine_quadratic_eqn<1>(__dr, __r, __ds, __s, __w);
	    else if (__w[0].index() == 2 && __w[1].index() == 2)
	      _S_refine_quadratic_eqn<2>(__dr, __r, __ds, __s, __w);
	  }
      }
    for (int __i = 0; __i < this->_M_order - 1; ++__i)
      this->_M_coeff[__i] = this->_M_b[__i];

    std::array<solution_t<_Real>, 2> __z2 = __quadratic(__s, __r, _Real{1});

    this->_M_add_zero(__z2[0]);
    this->_M_add_zero(__z2[1]);
  }

template<typename _Real>
  std::vector<solution_t<_Real>>
  _BairstowSolver<_Real>::solve()
  {
    this->_M_eps = _S_eps;

    this->_M_b[0] = _Real{1};
    this->_M_c[0] = _Real{1};

    while (this->_M_order > 2)
      this->_M_iterate();
    if (this->_M_order == 1)
      this->_M_add_zero(-this->_M_coeff[1]);

    return this->_M_zero;
  }

/**
 * Return the equations fouble by Bairstow's method.
 * There will be order/2 quadratic equations of the form
 * a[k] + a[k+1]x + x^2 = 0
 * and, if there is an odd term, a linear equation
 * a[order] + t = 0.
 */
template<typename _Real>
  std::vector<_Real>
  _BairstowSolver<_Real>::equations() const
  {
    std::vector<_Real> __eqs(this->_M_coeff.rbegin(), this->_M_coeff.rend());
    __eqs.erase(__eqs.begin() + __eqs.size() - 1, __eqs.end());
    __eqs[__eqs.size() - 1] *= _Real{-1};
    return __eqs;
  }

template<typename _Real>
  template<int _Index>
    void
    _BairstowSolver<_Real>::_S_refine_quadratic_eqn(_Real& __dr, _Real& __r,
					_Real& __ds, _Real& __s,
					std::array<solution_t<_Real>, 2>& __w)
    {
      const auto __p = std::experimental::make_array(__s, __r, _Real{1});
      __w[0] = __refine_solution_newton<3>(std::get<_Index>(__w[0]), __p);
      __w[1] = __refine_solution_newton<3>(std::get<_Index>(__w[1]), __p);
      const auto __sp = std::get<_Index>(__w[0]) * std::get<_Index>(__w[1]);
      const auto __rp = -(std::get<_Index>(__w[0]) + std::get<_Index>(__w[1]));
      __dr = std::real(__rp - __r);
      __r = std::real(__rp);
      __ds = std::real(__sp - __s);
      __s = std::real(__sp);
    }

} // namespace __gnu_cxx

template<typename _Real>
  void
  test_bairstow()
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);

    int order = 0;
    int MAX_TERMS = 100;
    while ((order < 2) || (order > MAX_TERMS - 1))
      {
	std::cout << "Polynomial order (2 - " << MAX_TERMS - 1 << "): ";
	std::cin >> order;
      }
    std::vector<_Real> a(order + 1);

    std::cout << "Enter coefficients, high order to low order.\n";
    for (int i = 0; i <= order; ++i)
      {
	std::cout << "a[" << i << "] = ";
	std::cin >> a[i];
      }

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

    std::cout << "\nSolution tests:\n";
    __gnu_cxx::_Polynomial<_Real> poly(a.begin(), a.end());
    for (const auto& z : zeros)
      {
	const auto idx = z.index();
	std::cout << "P(" << z << ") = ";
	if (idx == 1)
	  std::cout << poly(std::get<1>(z));
	else if (idx == 2)
	  std::cout << poly(std::get<2>(z));
	std::cout << '\n';
      }
  }

int
main()
{
  std::cout << "\ndouble\n======\n";
  test_bairstow<double>();

  std::cout << "\nlong double\n===========\n";
  test_bairstow<long double>();

  std::cout << "\nfloat\n=====\n";
  test_bairstow<float>();
}
