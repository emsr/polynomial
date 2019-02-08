// Math extensions -*- C++ -*-

// Copyright (C) 2018 Free Software Foundation, Inc.
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
 * @file static_polynomial.tcc Class outline definitions for static polynomial.
 */

/**
 * @def  _EXT_STATIC_POLYNOMIAL_TCC
 *
 * @brief  A guard for the static_polynomial class header.
 */
#ifndef _EXT_STATIC_POLYNOMIAL_TCC
#define _EXT_STATIC_POLYNOMIAL_TCC 1

#pragma GCC system_header

#if __cplusplus < 201402L
# include <bits/c++0x_warning.h>
#else

#include <iostream>

namespace __gnu_cxx //_GLIBCXX_VISIBILITY(default)
{

  /**
   * Write a polynomial to a stream.
   * The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp, std::size_t _Num>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& __os,
	       const _StaticPolynomial<_Tp, _Num>& __poly)
    {
      int __old_prec = __os.precision(std::numeric_limits<_Tp>::max_digits10);
      __os << "(";
      for (size_t __i = 0; __i < __poly.degree(); ++__i)
        __os << __poly.coefficient(__i) << ",";
      __os << __poly.coefficient(__poly.degree());
      __os << ")";
      __os.precision(__old_prec);
      return __os;
    }

} // namespace __gnu_cxx

#endif // C++14

#endif // _EXT_STATIC_POLYNOMIAL_TCC
