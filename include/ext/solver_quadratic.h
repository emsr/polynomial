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
 * @file solver_quadratic.h Class declaration for the quadratic solver.
 */

/**
 * @def  SOLVER_QUADRATIC_H
 *
 * @brief  A guard for the _QuadraticSolver class header.
 */
#ifndef SOLVER_QUADRATIC_H
#define SOLVER_QUADRATIC_H 1

#include <complex>

#include <ext/polynomial.h>
//#include <ext/solution.h> // For solution_t

namespace __gnu_cxx
{

  /**
   * A solver for complex-coefficient polynomials due to Laguerre.
   */
  template<typename _Real>
    class _QuadraticSolver
    {
    public:

      _QuadraticSolver(_Polynomial<std::complex<_Real>>& _P)
      : _M_poly(_P), _M_num_iters{0}
      { }

      //std::vector<solution_t<_Real>> solve();
      std::vector<std::complex<_Real>> solve();

      _Polynomial<std::complex<_Real>>
      step()
      {
	const auto __q = this->_M_root_quadratic();
	this->_M_poly.deflate(__q, _Real{10} * _S_eps);
	return __q;
      }

      int
      num_iters() const
      { return this->_M_num_iters; }

      int
      max_num_iters() const
      { return this->_M_max_num_iters; }

      const _Polynomial<std::complex<_Real>>&
      polynomial() const
      { return this->_M_poly; }

    private:

      // Estimated fractional roundoff error.
      static constexpr _Real _S_eps = std::numeric_limits<_Real>::epsilon();
      static constexpr _Real _S_tiny = _Real{10} * _S_eps;

      // Fractional roundoff error.
      _Real _M_eps = _Real{100} * _S_eps;

      int _M_max_iter = 50;

      _Polynomial<std::complex<_Real>> _M_root_quadratic();

      _Polynomial<std::complex<_Real>> _M_poly;

      int _M_num_iters = 0;
    };


} // namespace __gnu_cxx

#include <ext/solver_quadratic.tcc>

#endif // SOLVER_QUADRATIC_H
