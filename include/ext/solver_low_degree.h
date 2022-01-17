// Math extensions -*- C++ -*-

// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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
 * @file solver_low_degree.h
 *
 * This file is a GNU extension to the Standard C++ Library.
 *
 * This file contains the declarations of free functions for solving
 * quadratic, cubic, and quartic equations with real coefficients.
 */

/**
 * @def  SOLVER_LOW_DEGREE_H
 *
 * @brief  A guard for the low-degree polynomial solver functions header.
 */
#ifndef SOLVER_LOW_DEGREE_H
#define SOLVER_LOW_DEGREE_H 1

#include <array>

#include <ext/solution.h>

namespace emsr
{

  template<std::size_t Dim, typename Iter, typename NumTp>
    NumTp
    refine_solution_newton(NumTp z, const Iter& CC);

  template<std::size_t Dim, typename Iter, typename NumTp>
    NumTp
    refine_solution_halley(NumTp z, const Iter& CC);

  template<typename Real, typename Iter>
    std::array<solution_t<Real>, 2>
    quadratic(const Iter& coef);

  template<typename Real>
    inline std::array<solution_t<Real>, 2>
    quadratic(Real c0, Real c1, Real c2)
    {
      return quadratic<Real>(std::array<Real, 3>{c0, c1, c2});
    }

  template<typename Real, typename Iter>
    std::array<solution_t<Real>, 3>
    cubic(const Iter& coef);

  template<typename Real>
    inline std::array<solution_t<Real>, 3>
    cubic(Real c0, Real c1, Real c2, Real c3)
    {
      return cubic<Real>(std::array<Real, 4>{c0, c1, c2, c3});
    }

  template<typename Real, typename Iter>
    std::array<solution_t<Real>, 4>
    quartic(const Iter& coef);

  template<typename Real>
    inline std::array<solution_t<Real>, 4>
    quartic(Real c0, Real c1, Real c2, Real c3, Real c4)
    {
      return quartic<Real>(std::array<Real, 5>{c0, c1, c2, c3, c4});
    }

} // namespace emsr

#include <ext/solver_low_degree.tcc>

#endif // SOLVER_LOW_DEGREE_H
