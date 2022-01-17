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

namespace emsr
{

  /**
   * A solver for complex-coefficient polynomials due to Laguerre.
   */
  template<typename Real>
    class QuadraticSolver
    {
    public:

      QuadraticSolver(Polynomial<std::complex<Real>>& P)
      : m_poly(P), m_num_iters{0}
      { }

      //std::vector<solution_t<Real>> solve();
      std::vector<std::complex<Real>> solve();

      Polynomial<std::complex<Real>>
      step()
      {
	const auto q = this->m_root_quadratic();
	this->m_poly.deflate(q, Real{10} * s_eps);
	return q;
      }

      int
      num_iters() const
      { return this->m_num_iters; }

      int
      max_num_iters() const
      { return this->m_max_num_iters; }

      const Polynomial<std::complex<Real>>&
      polynomial() const
      { return this->m_poly; }

    private:

      // Estimated fractional roundoff error.
      static constexpr Real s_eps = std::numeric_limits<Real>::epsilon();
      static constexpr Real s_tiny = Real{10} * s_eps;

      // Fractional roundoff error.
      Real m_eps = Real{100} * s_eps;

      int m_max_iter = 50;

      Polynomial<std::complex<Real>> m_root_quadratic();

      Polynomial<std::complex<Real>> m_poly;

      int m_num_iters = 0;
    };

} // namespace emsr

#include <ext/solver_quadratic.tcc>

#endif // SOLVER_QUADRATIC_H
