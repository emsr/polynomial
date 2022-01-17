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

namespace emsr
{

  /**
   * Find a root of a complex-coefficient polynomial by Laguerre's method.
   * This routine can be iterated by dividing the original polynomial
   * by the root factor (z - x) where x is the found root and finding
   * the next root.
   */
  template<typename Real>
    std::complex<Real>
    LaguerreSolver<Real>::m_root_laguerre()
    {
      using cmplx = std::complex<Real>;

      const auto m = this->m_poly.degree();
      const int max_iter = this->m_max_iter();

      this->m_num_iters = 0;

      //if (m == 1)
	//return -this->m_poly.cefficient(0) / this->m_poly.cefficient(1);

      cmplx x{};
      for (int iter = 1; iter <= max_iter; ++iter)
	{
	  ++this->m_num_iters;

	  // Efficient computation of the polynomial
	  // and its first two derivatives. F stores P''(x)/2.
	  auto b = this->m_poly[m];
	  auto err = std::abs(b);
	  const auto abx = std::abs(x);
	  cmplx d{}, f{};
	  for (int j = m - 1; j >= 0; --j)
	    {
	      f = x * f + d;
	      d = x * d + b;
	      b = x * b + this->m_poly[j];
	      err = abx * err + std::abs(b);
	    }
	  err *= s_eps;
	  // Estimate of roundoff error in evaluating polynomial.
	  if (std::abs(b) <= err) // We have the root.
	    return x;

	  // Use Laguerre's formula.
	  const auto g = d / b;
	  const auto g2 = g * g;
	  const auto h = g2 - Real{2} * f / b;
	  const auto sq = std::sqrt(Real(m - 1)
				   * (Real(m) * h - g2));
	  auto gp = g + sq;
	  const auto gm = g - sq;
	  const auto abp = std::abs(gp);
	  const auto abm = std::abs(gm);
	  if (abp < abm)
	    gp = gm;
	  const auto dx = std::max(abp, abm) > Real{0}
			  ? Real(m) / gp
			  : std::polar(Real{1} + abx, Real(iter));
	  const auto x1 = x - dx;
	  if (x == x1)
	    return x;
	  if (iter % this->m_steps_per_frac != 0)
	    x = x1;
	  else
	    x -= s_frac[iter / this->m_steps_per_frac] * dx;
	}

      throw std::runtime_error("m_root_laguerre: Maximum number of iterations exceeded");
    }

  /**
   * Find all solutions of the  polynomial by stepping and 
   */
  template<typename Real>
    std::vector<std::complex<Real>>
    LaguerreSolver<Real>::solve()
    {
      using cmplx = std::complex<Real>;

      std::vector<cmplx> roots;
      const auto deg = this->m_poly.degree();
      roots.reserve(deg);
      for (unsigned i = 0; i < deg; ++i)
	{
	  const auto z0 = this->m_root_laguerre();

	  Polynomial<cmplx> zpoly({-z0, cmplx{1}});
	  this->m_poly /= zpoly;

	  roots.push_back(z0);
	}
      return roots;
    }

} // namespace emsr

#endif // SOLVER_LAGUERRE_TCC
