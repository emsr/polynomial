// Math extensions -*- C++ -*-

// Copyright (C) 2018-2019 Free Software Foundation, Inc.
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
 * @file solver_bairstow.h Class declaration for the Bairstow solver.
 */

/**
 * @def  SOLVER_BAIRSTOW_H
 *
 * @brief  A guard for the _BairstowSolver class header.
 */
#ifndef SOLVER_BAIRSTOW_H
#define SOLVER_BAIRSTOW_H 1

#include <ext/solver_low_degree.h>

namespace emsr
{

///
/// @todo: If you *don't* reverse the input array you solve for 1/z.
/// @todo Take a max_error.  If this->m_eps grows larger than max_error throw.
///
template<typename Real>
  class BairstowSolver
  {
  public:

    BairstowSolver(const std::vector<Real>& coeff,
		    unsigned int seed = std::random_device()())
    : m_coeff{coeff.rbegin(), coeff.rend()},
      m_b(coeff.size()), m_c(coeff.size()),
      m_order(coeff.size() - 1),
      m_urng(seed), m_pdf(Real{0}, Real{2})
    {
      if (this->m_coeff.size() == 0)
	throw std::domain_error("BairstowSolver: Coefficient size must be nonzero.");

      if (this->m_coeff[0] == Real{0})
	throw std::domain_error("BairstowSolver: Leading-order coefficient must be nonzero.");

      const auto scale = this->m_coeff[0];
      for (int i = 0; i <= this->m_order; ++i)
	this->m_coeff[i] /= scale;

      this->m_zero.reserve(this->m_coeff.size());
    }

    std::vector<solution_t<Real>> solve();
    std::vector<Real> equations() const;

  private:

    void m_iterate();

    template<int Index>
    static void s_refine_quadratic_eqn(Real& dr, Real& r,
					Real& ds, Real& s,
					std::array<solution_t<Real>, 2>& w);

    void
    m_add_zero(solution_t<Real> z)
    {
      this->m_zero.push_back(z);
      --this->m_order;
    }

    static constexpr Real s_ratio
	= Real{std::numeric_limits<Real>::digits}
	/ std::numeric_limits<double>::digits;
    static constexpr int s_max_rand_iter = 200 * s_ratio;
    static constexpr int s_max_error_iter = 500 * s_ratio;
    static constexpr auto s_eps_factor = Real{10} * s_ratio;
    static constexpr auto s_eps
      = s_eps_factor * std::numeric_limits<Real>::epsilon();

    std::vector<Real> m_coeff;
    std::vector<Real> m_b;
    std::vector<Real> m_c;
    std::vector<solution_t<Real>> m_zero;
    Real m_eps = s_eps;
    int m_order;
    bool m_precision_error = false;
    std::mt19937 m_urng;
    std::uniform_real_distribution<Real> m_pdf;
  };

} // namespace emsr

#include <ext/solver_bairstow.tcc>

#endif // SOLVER_BAIRSTOW_H
