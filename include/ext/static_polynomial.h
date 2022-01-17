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
 * @file static_polynomial.h Class definition for a static polynomial.
 */

/**
 * @def  STATIC_POLYNOMIAL_H
 *
 * @brief  A guard for the static_polynomial class header.
 */
#ifndef STATIC_POLYNOMIAL_H
#define STATIC_POLYNOMIAL_H 1

#include <limits>
#include <array>
#include <complex>
#include <iosfwd>

namespace emsr
{

  /**
   * This is a constant size polynomial.
   * It is really meant to just evaluate canned polynomial literals.
   */
  template<typename Tp, std::size_t Size>
    class StaticPolynomial
    {
    public:
      /**
       *  Typedefs.
       */
      using value_type = typename std::array<Tp, Size>::value_type;
      using reference = typename std::array<Tp, Size>::reference;
      using const_reference = typename std::array<Tp, Size>::const_reference;
      using pointer = typename std::array<Tp, Size>::pointer;
      using const_pointer = typename std::array<Tp, Size>::const_pointer;
      using iterator = typename std::array<value_type, Size>::iterator;
      using const_iterator = typename std::array<value_type, Size>::const_iterator;
      using reverse_iterator = typename std::array<value_type, Size>::reverse_iterator;
      using const_reverse_iterator = typename std::array<value_type, Size>::const_reverse_iterator;
      using size_type = typename std::array<Tp, Size>::size_type;
      using difference_type = typename std::array<Tp, Size>::difference_type;

      /**
       *  Create a zero degree polynomial with value zero.
       */
      constexpr
      StaticPolynomial()
      : m_coeff{}
      { }

      /**
       *  Copy ctor.
       */
      constexpr StaticPolynomial(const StaticPolynomial&) = default;
      constexpr StaticPolynomial(StaticPolynomial&&) = default;

      template<typename Up>
	constexpr
	StaticPolynomial(const StaticPolynomial<Up, Size>& poly)
	: m_coeff{}
	{
          for (auto i = 0ULL; i < Size; ++i)
	    this->m_coeff[i] = static_cast<value_type>(poly.m_coeff[i]);
	}

      /**
       *  Constructor from C-type array.
       */
      template<typename Up>
	constexpr
	StaticPolynomial(const Up (&arr)[Size])
	: m_coeff{}
	{
          for (auto i = 0ULL; i < Size; ++i)
	    this->m_coeff[i] = static_cast<value_type>(arr[i]);
	}

      /**
       *  Constructor from initializer_list array.
       */
      constexpr
      StaticPolynomial(std::initializer_list<Tp> il)
      : m_coeff{}
      {
	//static_assert(il.size() == Size, "");
	std::size_t i = 0;
	for (auto&& coeff : il)
	  this->m_coeff[i++] = coeff;
      }

      /**
       *  Create a polynomial - actually a monomial - of just one term.
       */
      constexpr explicit
      StaticPolynomial(value_type a, size_type degree = 0)
      : m_coeff(degree + 1)
      {
        static_assert(degree < Size, "StaticPolynomial: degree out of range");
        this->m_coeff[degree] = a;
      }

      /**
       *  Create a polynomial from an input iterator range of coefficients.
       */
      template<typename InIter,
	       typename = std::_RequireInputIter<InIter>>
	constexpr
	StaticPolynomial(const InIter& abegin, const InIter& aend)
	: m_coeff(abegin, aend)
	{ }

      /**
       *  Swap the polynomial with another polynomial.
       */
      void
      swap(StaticPolynomial& poly)
      { this->m_coeff.swap(poly.m_coeff); }

      /**
       *  Evaluate the polynomial at the input point.
       */
      constexpr value_type
      operator()(value_type x) const
      {
	if (this->degree() > 0)
	  {
	    value_type poly(this->coefficient(this->degree()));
	    for (int i = this->degree() - 1; i >= 0; --i)
	      poly = poly * x + this->coefficient(i);
	    return poly;
	  }
	else
	  return value_type{};
      }

      /**
       *  Evaluate the polynomial at the input point.
       */
      template<typename Tp2>
	constexpr auto
	operator()(Tp2 x) const
	-> decltype(value_type{} * Tp2())
	{
	  if (this->degree() > 0)
	    {
	      auto poly(this->coefficient(this->degree()) * Tp2(1));
	      for (int i = this->degree() - 1; i >= 0; --i)
		poly = poly * x + this->coefficient(i);
	      return poly;
	    }
	  else
	    return value_type{};
	}

      /**
       *  The following polynomial evaluations are done using
       *  a modified of Horner's rule which exploits the fact that
       *  the polynomial coefficients are all real.
       *  The algorithm is discussed in detail in:
       *  Knuth, D. E., The Art of Computer Programming: Seminumerical
       *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       *
       *  If n is the degree of the polynomial, n - 3 multiplies are
       *  saved and 4 * n - 6 additions are saved.
       */
      template<typename Tp2>
	constexpr auto
	operator()(std::complex<Tp2> z) const
	-> decltype(value_type{} * std::complex<Tp2>{})
	{
	  const auto r = Tp{2} * std::real(z);
	  const auto s = std::norm(z);
	  auto aa = this->coefficient(this->degree());
	  auto bb = this->coefficient(this->degree() - 1);
	  for (int j = 1; j <= this->degree(); ++j)
	    {
	      auto cc  = s * aa;
	      aa = bb + r * aa;
	      bb = this->coefficient(this->degree() - j) - cc;
	    }
	  return aa * z + bb;
	};

      /**
       *  Evaluate the polynomial at a range of input points.
       *  The output is written to the output iterator which
       *  must be large enough to contain the results.
       *  The next available output iterator is returned.
       */
      template<typename InIter, typename OutIter,
	       typename = std::_RequireInputIter<InIter>>
	constexpr OutIter
	operator()(const InIter& xbegin, const InIter& xend,
        	   OutIter& pbegin) const
	{
	  for (; xbegin != xend; ++xbegin)
	    pbegin++ = (*this)(xbegin++);
	  return pbegin;
	}

      //  Could/should this be done by output iterator range?
      template<size_type N>
	constexpr void
	eval(value_type x, std::array<value_type, N>& arr)
	{
	  if (arr.size() > 0)
	    {
	      arr.fill(value_type{});
	      const size_type sz = m_coeff.size();
	      arr[0] = this->coefficient(sz - 1);
              for (int i = sz - 2; i >= 0; --i)
		{
		  int nn = std::min(arr.size() - 1, sz - 1 - i);
		  for (int j = nn; j >= 1; --j)
		    arr[j] = arr[j] * x + arr[j - 1];
		  arr[0] = arr[0] * x + this->coefficient(i);
		}
	      //  Now put in the factorials.
	      value_type fact = value_type(1);
	      for (size_t i = 2; i < arr.size(); ++i)
		{
		  fact *= value_type(i);
		  arr[i] *= fact;
		}
	    }
	}

      /**
       *  Evaluate the polynomial and its derivatives at the point x.
       *  The values are placed in the output range starting with the
       *  polynomial value and continuing through higher derivatives.
       */
      template<typename OutIter>
	constexpr void
	eval(value_type x, OutIter b, OutIter e)
	{
	  if(b != e)
	    {
	      std::fill(b, e, value_type{});
	      constexpr size_type sz = m_coeff.size();
	      *b = m_coeff[sz - 1];
              for (int i = sz - 2; i >= 0; --i)
		{
		  for (auto it = std::reverse_iterator<OutIter>(e);
			   it != std::reverse_iterator<OutIter>(b) - 1; ++it)
		    *it = *it * x + *(it + 1);
		  *b = *b * x + m_coeff[i];
		}
	      //  Now put in the factorials.
	      int i = 0;
	      value_type fact = value_type(++i);
	      for (auto it = b + 1; it != e; ++it)
		{
		  fact *= value_type(i);
		  *it *= fact;
		  ++i;
		}
	    }
	}

      /**
       *  Evaluate the even part of the polynomial at the input point.
       */
      constexpr value_type
      eval_even(value_type x) const
      {
	if (this->degree() > 0)
	  {
	    auto odd = this->degree() % 2;
	    value_type poly(this->coefficient(this->degree() - odd));
	    for (int i = this->degree() - odd - 2; i >= 0; i -= 2)
	      poly = poly * x * x + this->coefficient(i);
	    return poly;
	  }
	else
	  return value_type{};
      }

      /**
       *  Evaluate the odd part of the polynomial at the input point.
       */
      constexpr value_type
      eval_odd(value_type x) const
      {
	if (this->degree() > 0)
	  {
	    auto even = (this->degree() % 2 == 0 ? 1 : 0);
	    value_type poly(this->coefficient(this->degree() - even));
	    for (int i = this->degree() - even - 2; i >= 0; i -= 2)
	      poly = poly * x * x + this->coefficient(i);
	    return poly * x;
	  }
	else
	  return value_type{};
      }

      /**
       *  Evaluate the even part of the polynomial using a modification
       *  of Horner's rule which exploits the fact that the polynomial
       *  coefficients are all real.
       *
       *  The algorithm is discussed in detail in:
       *  Knuth, D. E., The Art of Computer Programming: Seminumerical
       *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       *
       *  If n is the degree of the polynomial,
       *  n - 3 multiplies and 4 * n - 6 additions are saved.
       */
      template<typename Tp2>
	constexpr auto
	eval_even(std::complex<Tp2> z) const
	-> decltype(value_type{} * std::complex<Tp2>{})
	{
	  if (this->degree() > 0)
	    {
	      const auto zz = z * z;
	      const auto r = Tp{2} * std::real(zz);
	      const auto s = std::norm(zz);
	      auto odd = this->degree() % 2;
	      size_type n = this->degree() - odd;
	      auto aa = this->coefficient(n);
	      auto bb = this->coefficient(n - 2);
	      for (size_type j = 4; j <= n; j += 2)
		bb = this->coefficient(n - j)
		     - s * std::exchange(aa, bb + r * aa);
	      return aa * zz + bb;
	    }
	  else
	    return decltype(value_type{} * std::complex<Tp2>{}){};
	};

      /**
       *  Evaluate the odd part of the polynomial using a modification
       *  of Horner's rule which exploits the fact that the polynomial
       *  coefficients are all real.
       *
       *  The algorithm is discussed in detail in:
       *  Knuth, D. E., The Art of Computer Programming: Seminumerical
       *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       *
       *  If n is the degree of the polynomial,
       *  n - 3 multiplies and 4 * n - 6 additions are saved.
       */
      template<typename Tp2>
	constexpr auto
	eval_odd(std::complex<Tp2> z) const
	-> decltype(value_type{} * std::complex<Tp2>{})
	{
	  if (this->degree() > 0)
	    {
	      const auto zz = z * z;
	      const auto r = Tp{2} * std::real(zz);
	      const auto s = std::norm(zz);
	      auto even = (this->degree() % 2 == 0 ? 1 : 0);
	      size_type n = this->degree() - even;
	      auto aa = this->coefficient(n);
	      auto bb = this->coefficient(n - 2);
	      for (size_type j = 4; j <= n; j += 2)
		bb = this->coefficient(n - j)
		     - s * std::exchange(aa, bb + r * aa);
	      return z * (aa * zz + bb);
	    }
	  else
	    return decltype(value_type{} * std::complex<Tp2>{}){};
	};

      /**
       *  Return the derivative of the polynomial.
       */
      constexpr StaticPolynomial<Tp, (Size > 1 ? Size - 1 : 1)>
      derivative() const
      {
	StaticPolynomial<Tp, (Size > 1 ? Size - 1 : 1)> res;
	for (size_type i = 1; i <= this->degree(); ++i)
	  res.coefficient(i - 1, i * m_coeff[i]);
	return res;
      }

      /**
       *  Return the integral of the polynomial with given integration constant.
       */
      constexpr StaticPolynomial<Tp, Size + 1>
      integral(value_type c = value_type{}) const
      {
	StaticPolynomial<Tp, Size + 1> res;
	res.coefficient(0, c);
	for (size_type i = 0; i <= this->degree(); ++i)
	  res.coefficient(i + 1, m_coeff[i] / value_type(i + 1));
	return res;
      }

      /**
       * Unary plus.
       */
      constexpr StaticPolynomial
      operator+() const noexcept
      { return *this; }

      /**
       * Unary minus.
       */
      constexpr StaticPolynomial
      operator-() const
      { return StaticPolynomial(*this) *= value_type(-1); }

      /**
       *  Copy assignment.
       */
      constexpr StaticPolynomial&
      operator=(const StaticPolynomial&) = default;

      template<typename Up>
	StaticPolynomial&
	operator=(const StaticPolynomial<Up, Size>& poly)
	{
	  if (&poly != this)
	    {
	      this->m_coeff.clear();
	      for (const auto c : poly)
		this->m_coeff = static_cast<value_type>(c);
	      return *this;
	    }
	}

      /**
       *  Assign from an initialiser list.
       */
      constexpr StaticPolynomial&
      operator=(std::initializer_list<value_type> ila)
      {
	for (size_type i = 0;
	     i <= std::min(this->degree(), ila.size()); ++i)
	  this->m_coeff[i] = ila[i];
	return *this;
      }

      /**
       * Add a scalar to the polynomial.
       */
      StaticPolynomial&
      operator+=(const value_type& x)
      {
	this->m_coeff[0] += static_cast<value_type>(x);
	return *this;
      }

      /**
       * Subtract a scalar from the polynomial.
       */
      StaticPolynomial&
      operator-=(const value_type& x)
      {
	this->m_coeff[0] -= static_cast<value_type>(x);
	return *this;
      }

      /**
       * Multiply the polynomial by a scalar.
       */
      StaticPolynomial&
      operator*=(const value_type& c)
      {
	for (size_type i = 0; i < this->m_coeff.size(); ++i)
	  this->m_coeff[i] *= static_cast<value_type>(c);
	return *this;
      }

      /**
       * Divide the polynomial by a scalar.
       */
      StaticPolynomial&
      operator/=(const value_type& c)
      {
	for (size_type i = 0; i < this->m_coeff.size(); ++i)
	  this->m_coeff[i] /= static_cast<value_type>(c);
	return *this;
      }

      /**
       *  Return the degree or the power of the largest coefficient.
       */
      constexpr size_type
      degree() const
      { return (this->m_coeff.size() > 0 ? this->m_coeff.size() - 1 : 0); }

      /**
       * Return the @c ith coefficient with range checking.
       */
      constexpr value_type
      coefficient(size_type i) const
      { return this->m_coeff.at(i); }

      /**
       * Set coefficient @c i to @c val with range checking.
       */
      constexpr void
      coefficient(size_type i, value_type val)
      { this->m_coeff.at(i) = val; }

      /**
       * Return coefficient @c i.
       */
      constexpr value_type
      operator[](size_type i) const
      { return this->m_coeff[i]; }

      /**
       * Return coefficient @c i as an assignable quantity.
       */
      reference
      operator[](size_type i)
      { return this->m_coeff[i]; }

      /**
       * Return a const vector of coefficients.
       */
      constexpr const std::array<value_type, Size>
      coefficients() const noexcept
      { return this->m_coeff; }

      /**
       * Return a vector of coefficients.
       */
      constexpr std::array<value_type, Size>
      coefficients() noexcept
      { return this->m_coeff; }

      /**
       * Return a @c const pointer to the coefficient sequence.
       */
      constexpr const value_type*
      data() const noexcept
      { return this->m_coeff.data(); }

      /**
       * Return a @c pointer to the coefficient sequence.
       */
      constexpr value_type*
      data() noexcept
      { return this->m_coeff.data(); }

      iterator
      begin()
      { return this->m_coeff.begin(); }

      iterator
      end()
      { return this->m_coeff.end(); }

      const_iterator
      begin() const
      { return this->m_coeff.begin(); }

      const_iterator
      end() const
      { return this->m_coeff.end(); }

      const_iterator
      cbegin() const
      { return this->m_coeff.cbegin(); }

      const_iterator
      cend() const
      { return this->m_coeff.cend(); }

      reverse_iterator
      rbegin()
      { return this->m_coeff.rbegin(); }

      reverse_iterator
      rend()
      { return this->m_coeff.rend(); }

      const_reverse_iterator
      rbegin() const
      { return this->m_coeff.rbegin(); }

      const_reverse_iterator
      rend() const
      { return this->m_coeff.rend(); }

      const_reverse_iterator
      crbegin() const
      { return this->m_coeff.crbegin(); }

      const_reverse_iterator
      crend() const
      { return this->m_coeff.crend(); }

      template<typename Tp1>
	friend bool
	operator==(const StaticPolynomial<Tp1, Size>& pa,
		   const StaticPolynomial<Tp1, Size>& pb);

    private:

      std::array<value_type, Size> m_coeff;
    };

  /**
   *  Return true if two polynomials are equal.
   */
  template<typename Tp, std::size_t SizeA, std::size_t SizeB>
    inline constexpr bool
    operator==(const StaticPolynomial<Tp, SizeA>&,
	       const StaticPolynomial<Tp, SizeB>&)
    { return false; }

  template<typename Tp, std::size_t Size>
    inline constexpr bool
    operator==(const StaticPolynomial<Tp, Size>& pa,
	       const StaticPolynomial<Tp, Size>& pb)
    { return pa.m_coeff == pb.m_coeff; }

  /**
   *  Return false if two polynomials are equal.
   */
  template<typename Tp, std::size_t SizeA, std::size_t SizeB>
    inline constexpr bool
    operator!=(const StaticPolynomial<Tp, SizeA>& pa,
	       const StaticPolynomial<Tp, SizeB>& pb)
    { return true; }

  /**
   *  Return false if two polynomials are equal.
   */
  template<typename Tp, std::size_t Size>
    inline constexpr bool
    operator!=(const StaticPolynomial<Tp, Size>& pa,
	       const StaticPolynomial<Tp, Size>& pb)
    { return !(pa == pb); }

  /**
   * Return the sum of a polynomial with a scalar.
   */
  template<typename Tp, std::size_t Size>
    inline constexpr StaticPolynomial<Tp, Size>
    operator+(const StaticPolynomial<Tp, Size>& poly, const Tp& x)
    { return StaticPolynomial<Tp, Size>(poly) += x; }

  template<typename Tp, std::size_t Size>
    inline StaticPolynomial<Tp, Size>
    operator+(const Tp& x, const StaticPolynomial<Tp, Size>& poly)
    { return StaticPolynomial<Tp, Size>(poly) += x; }

  /**
   * Return the difference of a polynomial with a scalar.
   */
  template<typename Tp, std::size_t Size>
    inline constexpr StaticPolynomial<Tp, Size>
    operator-(const StaticPolynomial<Tp, Size>& poly, const Tp& x)
    { return StaticPolynomial<Tp, Size>(poly) -= x; }

  template<typename Tp, std::size_t Size>
    inline StaticPolynomial<Tp, Size>
    operator-(const Tp& x, const StaticPolynomial<Tp, Size>& poly)
    { return -StaticPolynomial<Tp, Size>(poly) += x; }

  /**
   * Return the product of a polynomial with a scalar.
   */
  template<typename Tp, std::size_t Size>
    inline constexpr StaticPolynomial<Tp, Size>
    operator*(const StaticPolynomial<Tp, Size>& poly, const Tp& x)
    { return StaticPolynomial<Tp, Size>(poly) *= x; }

  template<typename Tp, std::size_t Size>
    inline StaticPolynomial<Tp, Size>
    operator*(const Tp& x, const StaticPolynomial<Tp, Size>& poly)
    { return StaticPolynomial<Tp, Size>(poly) *= x; }

  /**
   * Return the quotient of a polynomial with a scalar.
   */
  template<typename Tp, std::size_t Size>
    inline constexpr StaticPolynomial<Tp, Size>
    operator/(const StaticPolynomial<Tp, Size>& poly, const Tp& x)
    { return StaticPolynomial<Tp, Size>(poly) /= x; }

  /**
   * Write a polynomial to a stream.
   * The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename Tp, std::size_t Size>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os,
	       const StaticPolynomial<Tp, Size>& poly);

  /**
   *  Return the sum of two polynomials.
   */
  template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
    inline constexpr StaticPolynomial<Tp, std::max(SizeP, SizeQ)>
    operator+(const StaticPolynomial<Tp, SizeP>& P,
	      const StaticPolynomial<Tp, SizeQ>& Q)
    {
      if constexpr (SizeP >= SizeQ)
	{
	  StaticPolynomial<Tp, SizeP> R = P;
	  for (std::size_t i = 0; i < SizeQ; ++i)
	    R[i] += Q[i];
	  return R;
	}
      else
	return Q + P;
    }

  /**
   *  Return the difference of two polynomials.
   */
  template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
    inline constexpr StaticPolynomial<Tp, std::max(SizeP, SizeQ)>
    operator-(const StaticPolynomial<Tp, SizeP>& P,
	      const StaticPolynomial<Tp, SizeQ>& Q)
    { return P + -Q; }

  /**
   *  Return the product of two polynomials.
   */
  template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
    inline constexpr StaticPolynomial<Tp, SizeP + SizeQ - 1>
    operator*(StaticPolynomial<Tp, SizeP> P,
	      StaticPolynomial<Tp, SizeQ> Q)
    {
      StaticPolynomial<Tp, P.degree() + Q.degree() + 1> R;
      for (std::size_t i = 0; i <= P.degree(); ++i)
	for (std::size_t j = 0; j <= Q.degree(); ++j)
	  R[i + j] = P[i] * Q[j];
      return R;
    }

  /**
   * Return the product of two polynomials.
   */
  template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
    inline constexpr StaticPolynomial<Tp, SizeP + SizeQ - 1>
    operator*(StaticPolynomial<Tp, SizeP> P,
	      StaticPolynomial<Tp, SizeQ> Q);

  /**
   * Return type for divmod.
   */
  template<typename Tp, std::size_t SizeN, std::size_t SizeD>
    struct divmod_t
    {
      static constexpr std::size_t
      SizeQuo = (SizeD <= SizeN) ? SizeN - SizeD + 1 : 1;
      static constexpr std::size_t
      SizeRem = (SizeD > 1) ? SizeD - 1 : 1;

      StaticPolynomial<Tp, SizeQuo> quo;
      StaticPolynomial<Tp, SizeRem> rem;
    };

  /**
   * Divide two polynomials returning the quotient and remainder.
   */
  template<typename Tp, std::size_t SizeN, std::size_t SizeD>
    constexpr divmod_t<Tp, SizeN, SizeD>
    divmod(StaticPolynomial<Tp, SizeN> num,
	   StaticPolynomial<Tp, SizeD> den);

  /**
   * Return the quotient of two polynomials.
   */
  template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
    inline constexpr StaticPolynomial<Tp,
		     divmod_t<Tp, SizeP, SizeQ>::SizeQuo>
    operator/(StaticPolynomial<Tp, SizeP> P,
	      StaticPolynomial<Tp, SizeQ> Q)
    { return divmod(P, Q).quo; }

  /**
   * Return the remainder of two polynomials.
   */
  template<typename Tp, std::size_t SizeP, std::size_t SizeQ>
    inline constexpr StaticPolynomial<Tp,
		     divmod_t<Tp, SizeP, SizeQ>::SizeRem>
    operator%(StaticPolynomial<Tp, SizeP> P,
	      StaticPolynomial<Tp, SizeQ> Q)
    { return divmod(P, Q).rem; }

} // namespace emsr

#include <ext/static_polynomial.tcc>

#endif // STATIC_POLYNOMIAL_H
