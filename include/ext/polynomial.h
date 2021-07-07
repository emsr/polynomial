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
 * @file polynomial.h Class declaration for a dense univariate polynomial.
 *
 * This file is a GNU extension to the Standard C++ Library.
 *
 * This file contains the declaration of a dense-polynomial class.
 */

/**
 * @def  _EXT_POLYNOMIAL_H
 *
 * @brief  A guard for the polynomial class header.
 */
#ifndef _EXT_POLYNOMIAL_H
#define _EXT_POLYNOMIAL_H 1

#pragma GCC system_header

#if __cplusplus < 201402L
# include <bits/c++0x_warning.h>
#else

#include <initializer_list>
#include <vector>
#include <iosfwd>
#include <limits>
#include <array>
#include <utility> // For exchange.
#include <type_traits>
#include <complex>

/**
 * This class is a dense univariate polynomial.
 *
 * This polynomial has a size of at least one - degree 0.  The zero polynomial
 * has a_0 = 0.  There is no null polynomial. size == degree + 1
 *
 * How to handle division?
 *   operator/ returns the quotient of two polynomials discarding the remainder.
 *   operator% returns the remainder of two polynomials discarding the quotient.
 *   The divmod functions computes both the quotient and the remainder.
 *   N. B. I could add the remquo function?
 * *
 * Largest coefficient:
 *   This class does not enforce the coefficient of largest power to be nonzero.
 *   The user concerned about this should use the deflate method with a
 *   tolerance for the minimim magnitude coefficient.
 *   This library provides another deflate method that divides out a polynomial
 *   that is expected to factor the polynomial within a tolerance on the
 *   minimim magnitude coefficient.  If the remainder is non-zero then an
 *   exception is thrown.
 *
 * Mixed-type evaluation:
 *   It is common to have integer coefficient polynomials that you wish to
 *   evaluate at floating opint or complex values, or real-valued  coefficient
 *   polynomials evaluate at complex values.  For this reason we need evaluation
 *   function templates in this library.
 *
 * It would be promote_t<complex::value_type> for complex.
 *
 */
namespace __gnu_cxx //_GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION


  template<typename, typename = std::void_t<>>
    struct __has_imag_t
    : std::false_type
    { };

  template<typename _Tp>
    struct __has_imag_t<_Tp, std::void_t<decltype(std::declval<_Tp&>().imag())>>
    : std::true_type
    { };

  template<typename _Tp>
    constexpr auto __has_imag_v = __has_imag_t<_Tp>::value;


  template<typename, typename = std::void_t<>>
    struct __has_value_type_t
    : std::false_type
    { };

  template<typename _Tp>
    struct __has_value_type_t<_Tp, std::void_t<typename _Tp::value_type>>
    : std::true_type
    { };

  template<typename _Tp>
    constexpr auto __has_value_type_v = __has_value_type_t<_Tp>::value;


  template<typename _Tp>
    class _Polynomial;

  template<typename _Tp>
    struct __real_type
    { using type = _Tp; };

  template<typename _Tp>
    struct __real_type<std::complex<_Tp>>
    { using type = _Tp; };

  template<typename _Tp>
    struct __real_type<_Polynomial<_Tp>>;

  template<typename _Tp>
    using __real_type_t = typename __real_type<_Tp>::type;


  /**
   * @brief A dense polynomial class with a contiguous array of coefficients.
   * The coefficients are lowest-order first:
   * @f[
   *    P(x) = a_0 + a_1 x + ... + a_n x^n
   * @f]
   */
  template<typename _Tp>
    class _Polynomial
    {
    public:
      /**
       * Typedefs.
       */
      using value_type = typename std::vector<_Tp>::value_type;
      using reference = typename std::vector<value_type>::reference;
      using const_reference = typename std::vector<value_type>::const_reference;
      using pointer = typename std::vector<value_type>::pointer;
      using const_pointer = typename std::vector<value_type>::const_pointer;
      using iterator = typename std::vector<value_type>::iterator;
      using const_iterator = typename std::vector<value_type>::const_iterator;
      using reverse_iterator
		= typename std::vector<value_type>::reverse_iterator;
      using const_reverse_iterator
		= typename std::vector<value_type>::const_reverse_iterator;
      using size_type = typename std::vector<value_type>::size_type;
      using difference_type = typename std::vector<value_type>::difference_type;
      using real_type = __real_type_t<_Tp>;

      /**
       * Create a zero degree polynomial with coefficient value zero.
       */
      _Polynomial()
      : _M_coeff(1)
      { }

      /**
       * Copy ctor.
       */
      _Polynomial(const _Polynomial&) = default;

      /**
       * Move ctor.
       */
      _Polynomial(_Polynomial&&) noexcept = default;

      template<typename _Up>
	_Polynomial(const _Polynomial<_Up>& __poly)
	: _M_coeff{}
	{
          for (const auto __c : __poly)
	    this->_M_coeff.push_back(static_cast<value_type>(__c));
          this->_M_set_scale();
	}

      /**
       * Create a monomial.
       */
      explicit
      _Polynomial(value_type __a, size_type __degree = 0)
      : _M_coeff(__degree + 1)
      { this->_M_coeff[__degree] = __a; }

      /**
       * Create a polynomial from an initializer list of coefficients.
       */
      _Polynomial(std::initializer_list<value_type> __ila)
      : _M_coeff(__ila)
      { this->_M_set_scale(); }

      /**
       * Create a polynomial from an input iterator range of coefficients.
       */
      template<typename InIter,
	       typename = std::_RequireInputIter<InIter>>
	_Polynomial(const InIter& __abegin, const InIter& __aend)
	: _M_coeff(__abegin, __aend)
	{ this->_M_set_scale(); }

      /**
       * Use Lagrange interpolation to construct a polynomial passing through
       * the data points.  The degree will be one less than the number of points.
       */
      template<typename InIter,
	       typename = std::_RequireInputIter<InIter>>
	_Polynomial(const InIter& __xbegin, const InIter& __xend,
		    const InIter& __ybegin)
	: _M_coeff()
	{
	  std::vector<_Polynomial<value_type>> __numer;
	  std::vector<_Polynomial<value_type>> __denom;
	  for (auto __xi = __xbegin; __xi != __xend; ++__xi)
	    {
	      for (auto __xj = __xi + 1; __xj != __xend; ++__xj)
		__denom.push_back(value_type(*__xj) - value_type(*__xi));
	      __numer.push_back({-value_type(*__xi), value_type{1}});
	    }
          this->_M_set_scale();
	}

      /**
       * Create a polynomial from a generator and a maximum degree.
       */
      template<typename Gen>
	_Polynomial(Gen __gen, size_type __degree)
	: _M_coeff()
	{
	  this->_M_coeff.reserve(__degree + 1);
	  for (size_type __k = 0; __k <= __degree; ++__k)
	    this->_M_coeff.push_back(__gen(__k));
          this->_M_set_scale();
	}

      /**
       * Swap the polynomial with another polynomial.
       */
      void
      swap(_Polynomial& __poly) noexcept
      { this->_M_coeff.swap(__poly._M_coeff); }

      /**
       * Evaluate the polynomial at the input point.
       */
      value_type
      operator()(value_type __x) const
      {
	if (this->degree() > 0)
	  {
	    if (std::abs(__x) <= real_type{1})
	      {
		value_type __poly(this->coefficient(this->degree()));
		for (int __i = this->degree() - 1; __i >= 0; --__i)
		  __poly = __poly * __x + this->coefficient(__i);
		return __poly;
	      }
	    else
	      {
		const auto __rx = real_type{1} / __x;
		value_type __poly(this->coefficient(0));
		for (int __i = 1; __i <= this->degree(); ++__i)
		  __poly = __poly * __rx + this->coefficient(__i);
		for (int __i = 1; __i <= this->degree(); ++__i)
		  __poly *= __x;
		return __poly;
	      }
	  }
	else
	  return this->coefficient(0);
      }

      /**
       * Evaluate the polynomial at the input point.
       */
      template<typename _Up>
	auto
	operator()(_Up __x) const
	-> decltype(value_type{} * _Up{})
	{
	  if (this->degree() > 0)
	    {
	      if (std::abs(__x) <= real_type{1})
		{
		  auto __poly(_Up{1} * this->coefficient(this->degree()));
		  for (int __i = this->degree() - 1; __i >= 0; --__i)
		    __poly = __poly * __x + this->coefficient(__i);
		  return __poly;
		}
	      else
		{
		  const auto __rx = real_type{1} / __x;
		  auto __poly(_Up{1} * this->coefficient(0));
		  for (int __i = 1; __i <= this->degree(); ++__i)
		    __poly = __poly * __rx + this->coefficient(__i);
		  for (int __i = 1; __i <= this->degree(); ++__i)
		    __poly *= __x;
		  return __poly;
		}
	    }
	  else
	    return _Up{1} * this->coefficient(0);
	}

      /**
       * Evaluate the polynomial using a modification of Horner's rule which
       * exploits the fact that the polynomial coefficients are all real.
       *
       * The algorithm is discussed in detail in:
       * Knuth, D. E., The Art of Computer Programming: Seminumerical
       * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       *
       * If n is the degree of the polynomial,
       * n - 3 multiplies and 4 * n - 6 additions are saved.
       */
      template<typename _Up>
	auto
	operator()(const std::complex<_Up>& __z) const
	-> std::enable_if_t<!__has_imag_v<_Tp>,
			    std::complex<std::decay_t<
		decltype(typename _Polynomial<_Tp>::value_type{} * _Up{})>>>;

      /**
       * Evaluate the polynomial at a range of input points.
       * The output is written to the output iterator which
       * must be large enough to contain the results.
       * The next available output iterator is returned.
       */
      template<typename InIter, typename OutIter,
	       typename = std::_RequireInputIter<InIter>>
	OutIter
	operator()(const InIter& __xbegin, const InIter& __xend,
        	   OutIter& __pbegin) const
	{
	  for (; __xbegin != __xend; ++__xbegin)
	    __pbegin++ = (*this)(__xbegin++);
	  return __pbegin;
	}

      template<size_type N>
	void
	eval(value_type __x, std::array<value_type, N>& __arr);

      /**
       * Evaluate the polynomial and its derivatives at the point x.
       * The values are placed in the output range starting with the
       * polynomial value and continuing through higher derivatives.
       */
      template<typename OutIter>
	void
	eval(value_type __x, OutIter __b, OutIter __e);

      /**
       * Evaluate the even part of the polynomial at the input point.
       */
      value_type
      eval_even(value_type __x) const;

      /**
       * Evaluate the odd part of the polynomial at the input point.
       */
      value_type
      eval_odd(value_type __x) const;

      /**
       * Evaluate the even part of the polynomial using a modification
       * of Horner's rule which exploits the fact that the polynomial
       * coefficients are all real.
       *
       * The algorithm is discussed in detail in:
       * Knuth, D. E., The Art of Computer Programming: Seminumerical
       * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       *
       * If n is the degree of the polynomial,
       * n - 3 multiplies and 4 * n - 6 additions are saved.
       */
      template<typename _Up>
	auto
	eval_even(const std::complex<_Up>& __z) const
	-> std::enable_if_t<!__has_imag_v<_Tp>,
			    std::complex<std::decay_t<
		decltype(typename _Polynomial<_Tp>::value_type{} * _Up{})>>>;

      /**
       * Evaluate the odd part of the polynomial using a modification
       * of Horner's rule which exploits the fact that the polynomial
       * coefficients are all real.
       *
       * The algorithm is discussed in detail in:
       * Knuth, D. E., The Art of Computer Programming: Seminumerical
       * Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       *
       * If n is the degree of the polynomial,
       * n - 3 multiplies and 4 * n - 6 additions are saved.
       */
      template<typename _Up>
	auto
	eval_odd(const std::complex<_Up>& __z) const
	-> std::enable_if_t<!__has_imag_v<_Tp>,
			std::complex<std::decay_t<
		decltype(typename _Polynomial<_Tp>::value_type{} * _Up{})>>>;

      /**
       * Return the derivative polynomial.
       */
      _Polynomial
      derivative() const
      {
	_Polynomial __res(value_type{},
			  this->degree() > 0UL ? this->degree() - 1 : 0UL);
	for (size_type __n = this->degree(), __i = 1; __i <= __n; ++__i)
	  __res._M_coeff[__i - 1] = __i * this->_M_coeff[__i];
	return __res;
      }

      /**
       * Return the derivative of the polynomial at the given point.
       */
      template<typename _Up>
        decltype(_Up{} * value_type{})
        derivative(_Up c) const
        {
	  using res_t = decltype(_Up{} * value_type{});
	  const int n = this->degree();
	  res_t res = real_type(n) * this->_M_coeff[n];
          for (int i = n - 1; i > 0; --i)
	    res = c * res + real_type(i) * this->_M_coeff[i];
	  return res;
        }

      /**
       * Return the integral polynomial with given integration constant.
       */
      _Polynomial
      integral(value_type __c = value_type{}) const
      {
	_Polynomial __res(value_type{}, this->degree() + 1);
	__res._M_coeff[0] = __c;
	for (size_type __n = this->degree(), __i = 0; __i <= __n; ++__i)
	  __res._M_coeff[__i + 1] = this->_M_coeff[__i] / value_type(__i + 1);
	return __res;
      }

      /**
       * Return the integral of the polynomial with given integration limits.
       */
      template<typename _Up>
        decltype(_Up{} * value_type{})
	integral(_Up a, _Up b) const
        {
	  using res_t = decltype(_Up{} * value_type{});
	  const int n = this->degree();
	  const auto coeff = this->_M_coeff[n] / real_type(n + 1);
	  res_t resa = coeff * a;
	  res_t resb = coeff * b;
	  for (int i = n - 1; i >= 0; --i)
	    {
	      const auto coeff = this->_M_coeff[i] / real_type(i + 1);
	      resa += coeff;
	      resa *= a;
	      resb += coeff;
	      resb *= b;
	    }
	  return resb - resa;
        }

      /**
       * Unary plus.
       */
      _Polynomial
      operator+() const noexcept
      { return *this; }

      /**
       * Unary minus.
       */
      _Polynomial
      operator-() const
      { return _Polynomial(*this) *= value_type(-1); }

      /**
       * Assign from a scalar.
       * The result is a zero degree polynomial equal to the scalar.
       */
      _Polynomial&
      operator=(const value_type& __x)
      {
	this->_M_coeff = {__x};
	return *this;
      }

      /**
       * Copy assignment.
       */
      _Polynomial&
      operator=(const _Polynomial&) = default;

      template<typename _Up>
	_Polynomial&
	operator=(const _Polynomial<_Up>& __poly)
	{
	  if (&__poly != this)
	    {
	      this->_M_coeff.clear();
	      for (const auto __c : __poly)
		this->_M_coeff.push_back(static_cast<value_type>(__c));
	      return *this;
	    }
	}

      /**
       * Assign from an initialiser list.
       */
      _Polynomial&
      operator=(std::initializer_list<value_type> __ila)
      {
	this->_M_coeff = __ila;
	return *this;
      }

      /**
       * Add a scalar to the polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator+=(const _Up& __x)
	{
	  this->_M_coeff[0] += static_cast<value_type>(__x);
	  return *this;
	}

      /**
       * Subtract a scalar from the polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator-=(const _Up& __x)
	{
	  this->_M_coeff[0] -= static_cast<value_type>(__x);
	  return *this;
	}

      /**
       * Multiply the polynomial by a scalar.
       */
      template<typename _Up>
	_Polynomial&
	operator*=(const _Up& __c)
	{
	  for (size_type __i = 0; __i < this->_M_coeff.size(); ++__i)
	    this->_M_coeff[__i] *= static_cast<value_type>(__c);
	  return *this;
	}

      /**
       * Divide the polynomial by a scalar.
       */
      template<typename _Up>
	_Polynomial&
	operator/=(const _Up& __c)
	{
	  for (size_type __i = 0; __i < this->_M_coeff.size(); ++__i)
	    this->_M_coeff[__i] /= static_cast<value_type>(__c);
	  return *this;
	}

      /**
       * Take the modulus of the polynomial relative to a scalar.
       * The result is always a zero polunomial.
       */
      template<typename _Up>
	_Polynomial&
	operator%=(const _Up&)
	{
	  this->degree(0UL); // Resize.
	  this->_M_coeff[0] = value_type{};
	  return *this;
	}

      /**
       * Add another polynomial to the polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator+=(const _Polynomial<_Up>& __poly)
	{
	  this->degree(std::max(this->degree(), __poly.degree()));
	  for (size_type __n = __poly.degree(), __i = 0; __i <= __n; ++__i)
	    this->_M_coeff[__i] += static_cast<value_type>(__poly._M_coeff[__i]);
	  return *this;
	}

      /**
       * Subtract another polynomial from the polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator-=(const _Polynomial<_Up>& __poly)
	{
	  // Resize if necessary.
	  this->degree(std::max(this->degree(), __poly.degree()));
	  for (size_type __n = __poly.degree(), __i = 0; __i <= __n; ++__i)
	    this->_M_coeff[__i] -= static_cast<value_type>(__poly._M_coeff[__i]);
	  return *this;
	}

      /**
       * Multiply the polynomial by another polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator*=(const _Polynomial<_Up>& __poly);

      /**
       * Divide the polynomial by another polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator/=(const _Polynomial<_Up>& __poly)
	{
	  _Polynomial<value_type >__quo, __rem;
	  divmod(*this, __poly, __quo, __rem);
	  *this = __quo;
	  return *this;
	}

      /**
       * Take the modulus of (modulate?) the polynomial relative to another polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator%=(const _Polynomial<_Up>& __poly)
	{
	  _Polynomial<value_type >__quo, __rem;
	  divmod(*this, __poly, __quo, __rem);
	  *this = __rem;
	  return *this;
	}

      /**
       * Shift the polynomial using the Horner scheme.
       * Given our polynomial
       * @f[
       *   P(x) = a_0 + a_1 x + a_2 x^2 + ...
       * @f]
       * Obtain a new polynomial
       * @f[
       *   Q(z) = P(x + s) = a_0 + a_1 (x + s) + a_2 (x + s)^2 + ... = b_0 + b_1 x + b_2 x^2
       * @f]
       */
      void
      shift(value_type shift)
      {
	if (shift == value_type{})
	  return;
        const int n = this->degree();
	for (int j = 1; j <= n; ++j)
	  for (int i = 1; i <= n - j + 1; ++i)
	    this->_M_coeff[n - i] += shift * this->_M_coeff[n - i + 1];
      }

      /**
       * Return the degree or the power of the largest coefficient.
       */
      size_type
      degree() const noexcept
      { return this->_M_coeff.size() - 1; }

      /**
       * Set the degree or the power of the largest coefficient.
       */
      void
      degree(size_type __degree)
      { this->_M_coeff.resize(__degree + 1UL); }

      /**
       * Return the size of the coefficient sequence.
       */
      size_type
      size() const noexcept
      { return this->_M_coeff.size(); }

      /**
       * Return the @c ith coefficient with range checking.
       */
      value_type
      coefficient(size_type __i) const
      { return this->_M_coeff.at(__i); }

      /**
       * Set coefficient @c i to @c val with range checking.
       */
      void
      coefficient(size_type __i, value_type __val)
      { this->_M_coeff.at(__i) = __val; }

      /**
       * Return coefficient @c i.
       */
      value_type
      operator[](size_type __i) const noexcept
      { return this->_M_coeff[__i]; }

      /**
       * Return coefficient @c i as an assignable quantity.
       */
      reference
      operator[](size_type __i) noexcept
      { return this->_M_coeff[__i]; }

      /**
       * Return a const vector of coefficients.
       */
      const std::vector<value_type>
      coefficients() const noexcept
      { return this->_M_coeff; }

      /**
       * Return a vector of coefficients.
       */
      std::vector<value_type>
      coefficients() noexcept
      { return this->_M_coeff; }

      /**
       * Return a @c const pointer to the coefficient sequence.
       */
      const value_type*
      data() const noexcept
      { return this->_M_coeff.data(); }

      /**
       * Return a @c pointer to the coefficient sequence.
       */
      value_type*
      data() noexcept
      { return this->_M_coeff.data(); }

      /**
       * Return an iterator to the beginning of the coefficient sequence.
       */
      iterator
      begin() noexcept
      { return this->_M_coeff.begin(); }

      /**
       * Return an iterator to one past the end of the coefficient sequence.
       */
      iterator
      end() noexcept
      { return this->_M_coeff.end(); }

      /**
       * Return a @c const iterator the beginning
       * of the coefficient sequence.
       */
      const_iterator
      begin() const noexcept
      { return this->_M_coeff.begin(); }

      /**
       * Return a @c const iterator to one past the end
       * of the coefficient sequence.
       */
      const_iterator
      end() const noexcept
      { return this->_M_coeff.end(); }

      /**
       * Return a @c const iterator the beginning
       * of the coefficient sequence.
       */
      const_iterator
      cbegin() const noexcept
      { return this->_M_coeff.cbegin(); }

      /**
       * Return a @c const iterator to one past the end
       * of the coefficient sequence.
       */
      const_iterator
      cend() const noexcept
      { return this->_M_coeff.cend(); }

      reverse_iterator
      rbegin() noexcept
      { return this->_M_coeff.rbegin(); }

      reverse_iterator
      rend() noexcept
      { return this->_M_coeff.rend(); }

      const_reverse_iterator
      rbegin() const noexcept
      { return this->_M_coeff.rbegin(); }

      const_reverse_iterator
      rend() const noexcept
      { return this->_M_coeff.rend(); }

      const_reverse_iterator
      crbegin() const noexcept
      { return this->_M_coeff.crbegin(); }

      const_reverse_iterator
      crend() const noexcept
      { return this->_M_coeff.crend(); }

      template<typename CharT, typename Traits, typename _Tp1>
	friend std::basic_istream<CharT, Traits>&
	operator>>(std::basic_istream<CharT, Traits>&, _Polynomial<_Tp1>&);

      template<typename _Tp1>
	friend bool
	operator==(const _Polynomial<_Tp1>& __pa,
		   const _Polynomial<_Tp1>& __pb);

      /**
       * Remove zero max-order coefficients.
       */
      _Polynomial&
      deflate(real_type __max_abs_coef)
      {
	size_type __n = this->degree();
	for (size_type __i = this->degree(); __i > 0; --__i)
	  if (std::abs(this->_M_coeff[__i]) < __max_abs_coef)
	    --__n;
	  else
	    break;
	this->degree(__n);
        return *this;
      }

      /**
       * Divide the polynomial by an input polynomia and remove zero
       * max-order coefficients.
       */
      _Polynomial&
      deflate(const _Polynomial<value_type>& __poly,
	      real_type __max_abs_coef)
      {
	_Polynomial<value_type> __quo, __rem;
	divmod(*this, __poly, __quo, __rem);

	// Remainder should be null.
	size_type __n = __rem.degree();
	for (size_type __i = __rem.degree(); __i > 0; --__i)
	  if (std::abs(__rem[__i]) < __max_abs_coef)
	    --__n;
	  else
	    break;

	if (__n == 0)
	  *this = __quo.deflate(__max_abs_coef);
	else
	  throw std::runtime_error("deflate: ");

        return *this;
      }

    private:

      /// Return the scale.
      real_type
      _M_get_scale() const
      { return this->_M_scale; }

      real_type _M_scale = real_type{1};

      void _M_set_scale();

      std::vector<value_type> _M_coeff;
    };

  template<typename _Tp>
    struct __real_type<_Polynomial<_Tp>>
    { using type = typename _Polynomial<_Tp>::real_type; };

  /**
   * Return the scale for a polynomial.
   */
  template<typename _Tp>
    get_scale(const _Polynomial<_Tp>& __poly)
    { return __poly._M_get_scale(); }

  /**
   * Return the scale for a number.
   */
  template<typename _Tp>
    get_scale(const _Tp& __x)
    { return std::abs(__x); }

  /**
   * Return the sum of a polynomial with a scalar.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() + _Up())>
    operator+(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() + _Up())>(__poly) += __x; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() + _Up())>
    operator+(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() + _Up())>(__x) += __poly; }

  /**
   * Return the difference of a polynomial with a scalar.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() - _Up())>
    operator-(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() - _Up())>(__poly) -= __x; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() - _Up())>
    operator-(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() - _Up())>(__x) -= __poly; }

  /**
   * Return the product of a polynomial with a scalar.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() * _Up())>
    operator*(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() * _Up())>(__poly) *= __x; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() * _Up())>
    operator*(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() * _Up())>(__x) *= __poly; }

  /**
   * Return the quotient of a polynomial with a scalar.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator/(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() / _Up())>(__poly) /= __x; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator%(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() / _Up())>(__poly) %= __x; }

  /**
   * Return the sum of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() + _Up())>
    operator+(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() + _Up())>(__pa) += __pb; }

  /**
   * Return the difference of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() - _Up())>
    operator-(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() - _Up())>(__pa) -= __pb; }

  /**
   * Return the product of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() * _Up())>
    operator*(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() * _Up())>(__pa) *= __pb; }

  /**
   * Return the quotient of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator/(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() / _Up())>(__pa) /= __pb; }

  /**
   * Return the modulus or remainder of one polynomial relative to another one.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator%(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() / _Up())>(__pa) %= __pb; }

  /**
   * Return the quotient of a scalar and a polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator/(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() / _Up())>(__x) /= __poly; }

  /**
   * Return the modulus or remainder of a scalar divided by a polynomial.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator%(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() / _Up())>(__x) %= __poly; }

  /**
   * Divide two polynomials returning the quotient and remainder.
   */
  template<typename _Tp>
    void
    divmod(const _Polynomial<_Tp>& __num, const _Polynomial<_Tp>& __den,
           _Polynomial<_Tp>& __quo, _Polynomial<_Tp>& __rem);

  /**
   * Write a polynomial to a stream.
   * The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& __os,
	       const _Polynomial<_Tp>& __poly);

  /**
   * Read a polynomial from a stream.
   * The input format can be a plain scalar (zero degree polynomial)
   * or a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& __is,
	       _Polynomial<_Tp>& __poly);

  /**
   * Return true if two polynomials are equal.
   */
  template<typename _Tp>
    inline bool
    operator==(const _Polynomial<_Tp>& __pa, const _Polynomial<_Tp>& __pb)
    { return __pa._M_coeff == __pb._M_coeff; }

  /**
   * Return false if two polynomials are equal.
   */
  template<typename _Tp>
    inline bool
    operator!=(const _Polynomial<_Tp>& __pa, const _Polynomial<_Tp>& __pb)
    { return !(__pa == __pb); }

  /**
   * See _Polynomial::swap().
   */
  template<typename _Tp>
    inline void
    swap(_Polynomial<_Tp>& __pa, _Polynomial<_Tp>& __pb)
    noexcept(noexcept(__pa.swap(__pb)))
    { __pa.swap(__pb); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include <ext/polynomial.tcc>

#endif // C++14

#endif // _EXT_POLYNOMIAL_H

