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

/**
 * @file horner.h Class declaration for Horner polynomial evaluation.
 *
 * This file is a GNU extension to the Standard C++ Library.
 */

/**
 * @def  _EXT_HORNER_H
 *
 * @brief  A guard for the horner functions header.
 */
#ifndef _EXT_HORNER_H
#define _EXT_HORNER_H 1

#pragma GCC system_header

#include <type_traits>

namespace __gnu_cxx //_GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

/**
 * Perform compile-time evaluation of a constant zero-order polynomial.
 */
template<typename _ArgT, typename _Coef0>
  constexpr std::conditional_t<std::is_integral<_ArgT>::value, double, _ArgT>
  horner(_ArgT __x, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral<_ArgT>::value,
					double, _ArgT>;
    return __arg_t{__c0};
  }

/**
 * Perform compile-time evaluation of a constant polynomial.
 * The polynomial coefficients are lowest-order first.
 */
template<typename _ArgT, typename _Coef0, typename... _Coef>
  constexpr std::conditional_t<std::is_integral<_ArgT>::value, double, _ArgT>
  horner(_ArgT __x, _Coef0 __c0, _Coef... __c)
  {
    using __arg_t = std::conditional_t<std::is_integral<_ArgT>::value,
					double, _ArgT>;
    return __arg_t{__c0} + __x * horner(__x, __c...);
  }


/**
 * Perform compile-time evaluation of a constant zero-order polynomial.
 * The polynomial coefficients are highest-order first.
 */
template<typename _ArgT, typename _Coef0>
  constexpr std::conditional_t<std::is_integral<_ArgT>::value, double, _ArgT>
  horner_big_end(_ArgT, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral<_ArgT>::value,
					double, _ArgT>;
    return __arg_t{__c0};
  }

/**
 * Perform compile-time evaluation of a constant first-order polynomial.
 * The polynomial coefficients are highest-order first.
 */
template<typename _ArgT, typename _Coef1, typename _Coef0>
  constexpr std::conditional_t<std::is_integral<_ArgT>::value, double, _ArgT>
  horner_big_end(_ArgT __x, _Coef1 __c1, _Coef0 __c0)
  {
    using __arg_t = std::conditional_t<std::is_integral<_ArgT>::value,
					double, _ArgT>;
    return horner_big_end(__x, __x * __arg_t{__c1} + __arg_t{__c0});
  }

/**
 * Perform compile-time evaluation of a constant polynomial.
 * The polynomial coefficients are highest-order first.
 */
template<typename _ArgT, typename _CoefN, typename _CoefNm1, typename... _Coef>
  constexpr std::conditional_t<std::is_integral<_ArgT>::value, double, _ArgT>
  horner_big_end(_ArgT __x, _CoefN __cn, _CoefNm1 __cnm1, _Coef... __c)
  {
    using __arg_t = std::conditional_t<std::is_integral<_ArgT>::value,
					double, _ArgT>;
    return horner_big_end(__x, __x * __arg_t{__cn} + __arg_t{__cnm1}, __c...);
  }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // _EXT_HORNER_H
