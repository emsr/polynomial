// Math extensions -*- C++ -*-

// Copyright (C) 2019 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/**
 * @file solver_laguerre.tcc Class declaration for the Laguerre solver.
 */

/**
 * @def  SOLVER_LAGUERRE_TCC
 *
 * @brief  A guard for the _LaguerreSolver class header.
 */
#ifndef SOLVER_LAGUERRE_TCC
#define SOLVER_LAGUERRE_TCC 1

namespace __gnu_cxx
{

  /**
   * Find a root of a complex-coefficient polynomial by Laguerre's method.
   * This routine can be iterated by dividing the original polynomial
   * by the root factor (z - x) where x is the found root and finding
   * the next root.
   */
  template<typename _Real>
    std::complex<_Real>
    _LaguerreSolver<_Real>::_M_root_laguerre()
    {
      using __cmplx = std::complex<_Real>;

      const auto __m = this->_M_poly.degree();
      const int __max_iter = this->_M_max_iter();

      this->_M_num_iters = 0;

      //if (__m == 1)
	//return -this->_M_poly.cefficient(0) / this->_M_poly.cefficient(1);

      __cmplx __x{};
      for (int __iter = 1; __iter <= __max_iter; ++__iter)
	{
	  ++this->_M_num_iters;

	  // Efficient computation of the polynomial
	  // and its first two derivatives. F stores P''(x)/2.
	  auto __b = this->_M_poly[__m];
	  auto __err = std::abs(__b);
	  const auto __abx = std::abs(__x);
	  __cmplx __d{}, __f{};
	  for (int __j = __m - 1; __j >= 0; --__j)
	    {
	      __f = __x * __f + __d;
	      __d = __x * __d + __b;
	      __b = __x * __b + this->_M_poly[__j];
	      __err = __abx * __err + std::abs(__b);
	    }
	  __err *= _S_eps;
	  // Estimate of roundoff error in evaluating polynomial.
	  if (std::abs(__b) <= __err) // We have the root.
	    return __x;

	  // Use Laguerre's formula.
	  const auto __g = __d / __b;
	  const auto __g2 = __g * __g;
	  const auto __h = __g2 - _Real{2} * __f / __b;
	  const auto __sq = std::sqrt(_Real(__m - 1)
				   * (_Real(__m) * __h - __g2));
	  auto __gp = __g + __sq;
	  const auto __gm = __g - __sq;
	  const auto __abp = std::abs(__gp);
	  const auto __abm = std::abs(__gm);
	  if (__abp < __abm)
	    __gp = __gm;
	  const auto __dx = std::max(__abp, __abm) > _Real{0}
			  ? _Real(__m) / __gp
			  : std::polar(_Real{1} + __abx, _Real(__iter));
	  const auto __x1 = __x - __dx;
	  if (__x == __x1)
	    return __x;
	  if (__iter % this->_M_steps_per_frac != 0)
	    __x = __x1;
	  else
	    __x -= _S_frac[__iter / this->_M_steps_per_frac] * __dx;
	}

      std::__throw_runtime_error(__N("_M_root_laguerre: "
				     "Maximum number of iterations exceeded"));
    }

  /**
   * Find all solutions of the  polynomial by stepping and 
   */
  template<typename _Real>
    std::vector<std::complex<_Real>>
    _LaguerreSolver<_Real>::solve()
    {
      using __cmplx = std::complex<_Real>;

      std::vector<__cmplx> __roots;
      const auto __deg = this->_M_poly.degree();
      __roots.reserve(__deg);
      for (unsigned __i = 0; __i < __deg; ++__i)
	{
	  const auto __z0 = this->_M_root_laguerre();

	  _Polynomial<__cmplx> __zpoly({-__z0, __cmplx{1}});
	  this->_M_poly /= __zpoly;

	  __roots.push_back(__z0);
	}
      return __roots;
    }

} // namespace __gnu_cxx

#endif // SOLVER_LAGUERRE_TCC
