#ifndef SOLVER_BAIRSTOW_TCC
#define SOLVER_BAIRSTOW_TCC 1

#include "solver_low_degree.h"

namespace __gnu_cxx
{

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

#endif // SOLVER_BAIRSTOW_TCC
