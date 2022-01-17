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
 * @file solver_quadratic.tcc Class declaration for the quadratic solver.
 */

/**
 * @def  SOLVER_QUADRATIC_TCC
 *
 * @brief  A guard for the _QuadraticSolver class header.
 */
#ifndef SOLVER_QUADRATIC_TCC
#define SOLVER_QUADRATIC_TCC 1

namespace emsr
{

  /**
   * I think this is trying to factor out a quadratic
   * from a complex-coefficient polynomial.
   */
  template<typename Tp>
    Polynomial<std::complex<Tp>>
    QuadraticSolver<Tp>::m_root_quadratic()
    {
      using Cmplx = std::complex<Tp>;
      using Poly = Polynomial<Cmplx>;

      if (this->m_poly.degree() <= 2)
	return this->m_poly;

      this->m_num_iters = 0;

      Cmplx c, b;
      Poly q, qq, rem;
      for (int iter = 0; iter < this->m_max_iter; ++iter)
	{
	  ++this->m_num_iters;

	  Poly d({c, b, Cmplx{1}});

	  // First division: r, s.
	  divmod(this->m_poly, d, q, rem);
	  const auto s = rem[0];
	  const auto r = rem[1];

	  // Second division: partial r, s with respect to c.
	  divmod(q, d, qq, rem);
	  const auto sc = -rem[0];
	  const auto rc = -rem[1];
	  const auto sb = -c * rc;
	  const auto rb = -b * rc + sc;

	  // Solve 2x2 equation.
	  const auto dv = Tp{1} / (sb * rc - sc * rb);
	  const auto delb = ( r * sc - s * rc) * dv;
	  b += delb;
	  const auto delc = (-r * sb + s * rb) * dv;
	  c += delc;
	  if ((std::abs(delb) <= this->m_eps * std::abs(b)
	      || std::abs(b) < s_tiny)
           && (std::abs(delc) <= this->m_eps * std::abs(c)
	      || std::abs(c) < s_tiny))
	    return Poly({c, b, Cmplx{1}});
	}
      throw std::runtime_error("m_root_quadratic: Maximum number of iterations exceeded");
    }

} // namespace emsr

#endif // SOLVER_QUADRATIC_TCC
