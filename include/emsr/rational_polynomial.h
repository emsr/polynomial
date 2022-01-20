
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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
 * @file rational_polynomial.h
 *
 * This file contains the declaration of a ratio of two polynomials.
 * @see polynomial.h
 */

/**
 * @def  RATIONAL_POLYNOMIAL_H
 *
 * @brief  A guard for the rational_polynomial class header.
 */
#ifndef RATIONAL_POLYNOMIAL_H
#define RATIONAL_POLYNOMIAL_H 1

#include <iostream>

#include <emsr/polynomial.h>

namespace emsr
{

  /**
   *
   */
  template<typename Tp>
    class RationalPolynomial
    {
    public:
      /**
       * Typedefs.
       */
      using polynomial_type = Polynomial<Tp>;
      using value_type = typename polynomial_type::value_type;
// The above should be Polynomial<Tp>to be consistent with rational.h
/* These might not make sense.
      using reference = typename polynomial_type::reference;
      using const_reference = typename polynomial_type::const_reference;
      using pointer = typename polynomial_type::pointer;
      using const_pointer = typename polynomial_type::const_pointer;
      using iterator = typename polynomial_type::iterator;
      using const_iterator = typename polynomial_type::const_iterator;
      using reverse_iterator = typename polynomial_type::reverse_iterator;
      using const_reverse_iterator = typename polynomial_type::const_reverse_iterator;
*/
      using size_type = typename polynomial_type::size_type;
      using difference_type = typename polynomial_type::difference_type;

      /**
       * Create a zero degree polynomial with value zero.
       */
      RationalPolynomial()
      : m_num(), m_den()
      { }

      /**
       * Copy ctor.
       */
      RationalPolynomial(const RationalPolynomial&) = default;

      /**
       * Move ctor.
       */
      RationalPolynomial(RationalPolynomial&&) = default;

      RationalPolynomial(const Polynomial<Tp>& num,
      			  const Polynomial<Tp>& den)
      : m_num(num), m_den(den)
      { }

      explicit RationalPolynomial(const Polynomial<Tp>& num)
      : m_num(num), m_den(Polynomial<Tp>(Tp{1}))
      { }

      /**
       * Evaluate the polynomial at the input point.
       */
      value_type
      operator()(value_type x) const
      { return this->m_num(x) / this->m_den(x); }

      /**
       * Unary plus.
       */
      RationalPolynomial
      operator+() const
      { return *this; }

      /**
       * Unary minus.
       */
      RationalPolynomial
      operator-() const
      { return RationalPolynomial(*this) *= value_type(-1); }

      /**
       * Copy assignment.
       */
      RationalPolynomial&
      operator=(const RationalPolynomial&) = default;

      /**
       * Move assignment.
       */
      RationalPolynomial&
      operator=(RationalPolynomial&&) = default;

      /**
       * Add a rational polynomial to this rational polynomial.
       */
      RationalPolynomial&
      operator+=(const RationalPolynomial& x)
      {
        this->numer() = this->numer() * x.denom() + this->denom() * x.numer();
        this->denom() *= x.denom();
	return *this;
      }

      /**
       * Subtract a rational polynomial from this rational polynomial.
       */
      RationalPolynomial&
      operator-=(const RationalPolynomial& x)
      {
        this->numer() = this->numer() * x.denom() - this->denom() * x.numer();
        this->denom() *= x.denom();
	return *this;
      }

      /**
       * Multiply this rational polynomial by a rational polynomial.
       */
      RationalPolynomial&
      operator*=(const RationalPolynomial& x)
      {
	this->numer() *= x.numer();
	this->denom() *= x.denom();
	return *this;
      }

      /**
       * Divide this rational polynomial by a rational polynomial.
       */
      RationalPolynomial&
      operator/=(const RationalPolynomial& x)
      {
	this->numer() *= x.denom();
	this->denom() *= x.numer();
	return *this;
      }

      const Polynomial<value_type>&
      numer() const
      { return this->m_num; }

      Polynomial<value_type>&
      numer()
      { return this->m_num; }

      const Polynomial<value_type>&
      denom() const
      { return this->m_den; }

      Polynomial<value_type>&
      denom()
      { return this->m_den; }

    private:

      Polynomial<value_type> m_num;
      Polynomial<value_type> m_den;
    };

  /**
   * Write a polynomial to a stream.
   * The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename Tp>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os, const RationalPolynomial<Tp>& poly)
    {
      os << poly.numer() << "/" << poly.denom();
      return os;
    }

  /**
   * Read a polynomial from a stream.
   * The input format can be a plain scalar (zero degree polynomial)
   * or a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename Tp>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& is, RationalPolynomial<Tp>& poly)
    {
      Polynomial<Tp> numer, denom;
      is >> numer;
      if (!is.fail())
	{
	  CharT ch;
	  is >> ch;
	  if (ch != '/')
	    is.setstate(std::ios_base::failbit);
	  else
	    {
	      is >> denom;
	      if (!is.fail())
		{
		  poly.numer() = numer;
		  poly.denom() = denom;
		}
	    }
	}
      return is;
    }

} // namespace emsr

#endif // RATIONAL_POLYNOMIAL_H
