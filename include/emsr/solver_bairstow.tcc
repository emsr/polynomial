
// Copyright (C) 2018-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

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
 * @file solver_bairstow.tcc Class outline definitions for the Bairstow solver.
 */

/**
 * @def  SOLVER_BAIRSTOW_TCC
 *
 * @brief  A guard for the BairstowSolver class header.
 */
#ifndef SOLVER_BAIRSTOW_TCC
#define SOLVER_BAIRSTOW_TCC 1

#include <emsr/solver_low_degree.h>

namespace emsr
{

/**
 * Attempt to find a quadratic factor of the polynomial by Bairstow's method.
 */
template<typename Real>
  void
  BairstowSolver<Real>::m_iterate()
  {
    auto r = Real{0};
    auto s = Real{0};
    auto dr = Real{1};
    auto ds = Real{0};

    int iter = 1;
    this->m_precision_error = false;

    while (std::abs(dr) + std::abs(ds) > this->m_eps)
      {
	if (iter % s_max_rand_iter == 0)
	  r = this->m_pdf(this->m_urng);
	if (iter % s_max_error_iter == 0)
	  {
	    this->m_eps *= s_eps_factor;
	    this->m_precision_error = true;
	    std::cout << "Loss of precision: " << this->m_eps << '\n';
	  }

	this->m_b[1] = this->m_coeff[1] - r;
	this->m_c[1] = this->m_b[1] - r;
	for (int i = 2; i <= this->m_order; ++i)
	  {
	    this->m_b[i] = this->m_coeff[i]
			- r * this->m_b[i - 1] - s * this->m_b[i - 2];
	    this->m_c[i] = this->m_b[i]
			- r * this->m_c[i - 1] - s * this->m_c[i - 2];
	  }

	auto dn = this->m_c[this->m_order - 1]
		  * this->m_c[this->m_order - 3]
		  - this->m_c[this->m_order - 2]
		  * this->m_c[this->m_order - 2];
	if (std::abs(dn) < std::numeric_limits<Real>::epsilon())
	  {
	    dr = Real{1};
	    ds = Real{1};
	  }
	else
	  {
	    auto drn = this->m_b[this->m_order]
		       * this->m_c[this->m_order - 3]
		       - this->m_b[this->m_order - 1]
		       * this->m_c[this->m_order - 2];
	    auto dsn = this->m_b[this->m_order - 1]
		       * this->m_c[this->m_order - 1]
		       - this->m_b[this->m_order]
		       * this->m_c[this->m_order - 2];
	    dr = drn / dn;
	    ds = dsn / dn;
	  }

	r += dr;
	s += ds;
	++iter;
	if (std::abs(dr) + std::abs(ds) <= this->m_eps)
	  {
	    // Before exiting, solve the quadratic, refine the roots,
	    // multiply out the resulting polynomial, extract
	    // the new coefficients, get new dr, r, ds, s.
	    auto w = quadratic(s, r, Real{1});
	    if (w[0].index() == 1 && w[1].index() == 1)
	      s_refine_quadratic_eqn<1>(dr, r, ds, s, w);
	    else if (w[0].index() == 2 && w[1].index() == 2)
	      s_refine_quadratic_eqn<2>(dr, r, ds, s, w);
	  }
      }
    for (int i = 0; i < this->m_order - 1; ++i)
      this->m_coeff[i] = this->m_b[i];

    std::array<solution_t<Real>, 2> z2 = quadratic(s, r, Real{1});

    this->m_add_zero(z2[0]);
    this->m_add_zero(z2[1]);
  }

template<typename Real>
  std::vector<solution_t<Real>>
  BairstowSolver<Real>::solve()
  {
    this->m_eps = s_eps;

    this->m_b[0] = Real{1};
    this->m_c[0] = Real{1};

    while (this->m_order > 2)
      this->m_iterate();
    if (this->m_order == 1)
      this->m_add_zero(-this->m_coeff[1]);

    return this->m_zero;
  }

/**
 * Return the equations fouble by Bairstow's method.
 * There will be order/2 quadratic equations of the form
 * a[k] + a[k+1]x + x^2 = 0
 * and, if there is an odd term, a linear equation
 * a[order] + t = 0.
 */
template<typename Real>
  std::vector<Real>
  BairstowSolver<Real>::equations() const
  {
    std::vector<Real> eqs(this->m_coeff.rbegin(), this->m_coeff.rend());
    eqs.erase(eqs.begin() + eqs.size() - 1, eqs.end());
    eqs[eqs.size() - 1] *= Real{-1};
    return eqs;
  }

template<typename Real>
  template<int Index>
    void
    BairstowSolver<Real>::s_refine_quadratic_eqn(Real& dr, Real& r,
					Real& ds, Real& s,
					std::array<solution_t<Real>, 2>& w)
    {
      const auto p = std::array<Real, 3>{s, r, Real{1}};
      w[0] = refine_solution_newton<3>(std::get<Index>(w[0]), p);
      w[1] = refine_solution_newton<3>(std::get<Index>(w[1]), p);
      const auto sp = std::get<Index>(w[0]) * std::get<Index>(w[1]);
      const auto rp = -(std::get<Index>(w[0]) + std::get<Index>(w[1]));
      dr = std::real(rp - r);
      r = std::real(rp);
      ds = std::real(sp - s);
      s = std::real(sp);
    }

} // namespace emsr

#endif // SOLVER_BAIRSTOW_TCC
