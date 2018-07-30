// Math extensions -*- C++ -*-

// Copyright (C) 2016-2018 Free Software Foundation, Inc.
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

#ifndef SOLVER_LOW_DEGREE_H
#define SOLVER_LOW_DEGREE_H 1

/**
 * @file solver_low_degree.h
 *
 * This file is a GNU extension to the Standard C++ Library.
 *
 * This file contains the declarations of free functions for solving
 * quadratic, cubic, and quartic equations with real coefficients.
 */

#pragma GCC system_header

#if __cplusplus < 201703L
# include <bits/c++0x_warning.h>
#else

#include <experimental/array>

#include "solution.h"

namespace __gnu_cxx
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 2>
    __quadratic(const _Iter& __coef);

  template<typename _Real>
    inline std::array<solution_t<_Real>, 2>
    __quadratic(_Real __c0, _Real __c1, _Real __c2)
    {
      using std::experimental::make_array;
      return __quadratic<_Real>(make_array(__c0, __c1, __c2));
    }

  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 3>
    __cubic(const _Iter& __coef);

  template<typename _Real>
    inline std::array<solution_t<_Real>, 3>
    __cubic(_Real __c0, _Real __c1, _Real __c2, _Real __c3)
    {
      using std::experimental::make_array;
      return __cubic<_Real>(make_array(__c0, __c1, __c2, __c3));
    }

  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 4>
    __quartic(const _Iter& __coef);

  template<typename _Real>
    inline std::array<solution_t<_Real>, 4>
    __quartic(_Real __c0, _Real __c1, _Real __c2, _Real __c3, _Real __c4)
    {
      using std::experimental::make_array;
      return __quartic<_Real>(make_array(__c0, __c1, __c2, __c3, __c4));
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include "solver_low_degree.tcc"

#endif // C++17

#endif // SOLVER_LOW_DEGREE_H
