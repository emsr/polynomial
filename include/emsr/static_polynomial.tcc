
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
 * @file static_polynomial.tcc Class outline definitions for static polynomial.
 */

/**
 * @def  STATIC_POLYNOMIAL_TCC
 *
 * @brief  A guard for the static_polynomial class header.
 */
#ifndef STATIC_POLYNOMIAL_TCC
#define STATIC_POLYNOMIAL_TCC 1

#include <iostream>

namespace emsr
{

  /**
   * Write a polynomial to a stream.
   * The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename Tp, std::size_t Size>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os,
	       const StaticPolynomial<Tp, Size>& poly)
    {
      int old_prec = os.precision(std::numeric_limits<Tp>::max_digits10);
      os << "(";
      for (size_t i = 0; i < poly.degree(); ++i)
        os << poly.coefficient(i) << ",";
      os << poly.coefficient(poly.degree());
      os << ")";
      os.precision(old_prec);
      return os;
    }

  /**
   * Divide two polynomials returning the quotient and remainder.
   */
  template<typename Tp, std::size_t SizeN, std::size_t SizeD>
    constexpr divmod_t<Tp, SizeN, SizeD>
    divmod(StaticPolynomial<Tp, SizeN> num,
	   StaticPolynomial<Tp, SizeD> den)
    {
      constexpr auto DegN = num.degree();
      constexpr auto DegD = den.degree();
      auto rem = num;
      auto quo = StaticPolynomial<Tp, SizeN>{};
      if (DegD <= DegN)
	{
	  for (std::ptrdiff_t k = DegN - DegD; k >= 0; --k)
	    {
	      quo.coefficient(k, rem.coefficient(DegD + k)
				   / den.coefficient(DegD));
	      for (int j = DegD + k - 1; j >= k; --j)
		rem.coefficient(j, rem.coefficient(j)
				       - quo.coefficient(k)
				       * den.coefficient(j - k));
	    }
	}
      divmod_t<Tp, SizeN, SizeD> ret;
      for (std::size_t i = 0ULL; i < divmod_t<Tp, SizeN, SizeD>::SizeQuo; ++i)
        ret.quo[i] = quo[i];
      for (std::size_t i = 0ULL; i < divmod_t<Tp, SizeN, SizeD>::SizeRem; ++i)
        ret.rem[i] = rem[i];
      return ret;
    }

} // namespace emsr

#endif // STATIC_POLYNOMIAL_TCC
