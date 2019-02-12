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
 * @file solver_laguerre.h Class declaration for the Laguerre solver.
 */

/**
 * @def  SOLVER_LAGUERRE_H
 *
 * @brief  A guard for the _LaguerreSolver class header.
 */
#ifndef SOLVER_LAGUERRE_H
#define SOLVER_LAGUERRE_H 1

#include <complex>

#include <ext/polynomial.h>

namespace __gnu_cxx
{

  /**
   * A solver for complex-coefficient polynomials due to Laguerre.
   */
  template<typename _Real>
    class _LaguerreSolver
    {
    public:

      _LaguerreSolver(_Polynomial<std::complex<_Real>>& _P)
      : _M_poly(_P), _M_num_iters{0}
      { }

      std::vector<std::complex<_Real>> solve();

      std::complex<_Real>
      step()
      {
	using __cmplx = std::complex<_Real>;

	const auto __z0 = this->_M_root_laguerre();

	_Polynomial<__cmplx> __zpoly({-__z0, __cmplx{1}});
	this->_M_poly /= __zpoly;

	return __z0;
      }

      int
      num_iters() const
      { return this->_M_num_iters; }

      int
      max_num_iters() const
      { return this->_M_max_iter(); }

      int
      num_steps_per_frac() const
      { return this->_M_steps_per_frac; }

      const _Polynomial<std::complex<_Real>>&
      polynomial() const
      { return this->_M_poly; }

      _LaguerreSolver&
      num_steps_per_frac(int num)
      {
	this->_M_steps_per_frac = num;
	return *this;
      }

    private:

      // Estimated fractional roundoff error.
      static constexpr _Real _S_eps = std::numeric_limits<_Real>::epsilon();

      // Number of fractional values.
      static constexpr int _S_num_fracs = 8;
      // Fractions used to break a limit cycle (in a heap).
      static constexpr _Real
      _S_frac[_S_num_fracs + 1]
      {0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 1.0};

      // Number of steps taken before trying a new fraction.
      int _M_steps_per_frac = 10;

      int
      _M_max_iter()
      { return this->_M_steps_per_frac * _S_num_fracs; }

      std::complex<_Real> _M_root_laguerre();

      _Polynomial<std::complex<_Real>> _M_poly;

      int _M_num_iters = 0;
    };

} // namespace __gnu_cxx

#include <ext/solver_laguerre.tcc>

#endif // SOLVER_LAGUERRE_H
