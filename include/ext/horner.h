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
 * @file horner.h Class declaration for Horner polynomial evaluation.
 *
 * This file is a GNU extension to the Standard C++ Library.
 */

/**
 * @def  HORNER_H
 *
 * @brief  A guard for the horner functions header.
 */
#ifndef HORNER_H
#define HORNER_H 1

#pragma GCC system_header

#include <type_traits>

namespace emsr
{

/**
 * Perform compile-time evaluation of a constant zero-order polynomial.
 */
template<typename ArgT, typename Coef0>
  constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
  horner(ArgT x, Coef0 c0)
  {
    using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
					double, ArgT>;
    return arg_t{c0};
  }

/**
 * Perform compile-time evaluation of a constant polynomial.
 * The polynomial coefficients are lowest-order first.
 */
template<typename ArgT, typename Coef0, typename... Coef>
  constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
  horner(ArgT x, Coef0 c0, Coef... c)
  {
    using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
					double, ArgT>;
    return arg_t{c0} + x * horner(x, c...);
  }


/**
 * Perform compile-time evaluation of a constant zero-order polynomial.
 * The polynomial coefficients are highest-order first.
 */
template<typename ArgT, typename Coef0>
  constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
  horner_big_end(ArgT, Coef0 c0)
  {
    using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
					double, ArgT>;
    return arg_t{c0};
  }

/**
 * Perform compile-time evaluation of a constant first-order polynomial.
 * The polynomial coefficients are highest-order first.
 */
template<typename ArgT, typename Coef1, typename Coef0>
  constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
  horner_big_end(ArgT x, Coef1 c1, Coef0 c0)
  {
    using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
					double, ArgT>;
    return horner_big_end(x, x * arg_t{c1} + arg_t{c0});
  }

/**
 * Perform compile-time evaluation of a constant polynomial.
 * The polynomial coefficients are highest-order first.
 */
template<typename ArgT, typename CoefN, typename CoefNm1, typename... Coef>
  constexpr std::conditional_t<std::is_integral<ArgT>::value, double, ArgT>
  horner_big_end(ArgT x, CoefN cn, CoefNm1 cnm1, Coef... c)
  {
    using arg_t = std::conditional_t<std::is_integral<ArgT>::value,
					double, ArgT>;
    return horner_big_end(x, x * arg_t{cn} + arg_t{cnm1}, c...);
  }

} // namespace emsr

#endif // HORNER_H
