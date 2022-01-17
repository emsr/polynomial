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

namespace emsr
{

  /**
   * A solver for complex-coefficient polynomials due to Laguerre.
   */
  template<typename Real>
    class LaguerreSolver
    {
    public:

      LaguerreSolver(Polynomial<std::complex<Real>>& P)
      : m_poly(P), m_num_iters{0}
      { }

      std::vector<std::complex<Real>> solve();

      std::complex<Real>
      step()
      {
	using cmplx = std::complex<Real>;

	const auto z0 = this->m_root_laguerre();

	Polynomial<cmplx> zpoly({-z0, cmplx{1}});
	this->m_poly /= zpoly;

	return z0;
      }

      int
      num_iters() const
      { return this->m_num_iters; }

      int
      max_num_iters() const
      { return this->m_max_iter(); }

      int
      num_steps_per_frac() const
      { return this->m_steps_per_frac; }

      const Polynomial<std::complex<Real>>&
      polynomial() const
      { return this->m_poly; }

      LaguerreSolver&
      num_steps_per_frac(int num)
      {
	this->m_steps_per_frac = num;
	return *this;
      }

    private:

      // Estimated fractional roundoff error.
      static constexpr Real s_eps = std::numeric_limits<Real>::epsilon();

      // Number of fractional values.
      static constexpr int s_num_fracs = 8;
      // Fractions used to break a limit cycle (in a heap).
      static constexpr Real
      s_frac[s_num_fracs + 1]
      {0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 1.0};

      // Number of steps taken before trying a new fraction.
      int m_steps_per_frac = 10;

      int
      m_max_iter()
      { return this->m_steps_per_frac * s_num_fracs; }

      std::complex<Real> m_root_laguerre();

      Polynomial<std::complex<Real>> m_poly;

      int m_num_iters = 0;
    };

} // namespace emsr

#include <ext/solver_laguerre.tcc>

#endif // SOLVER_LAGUERRE_H
