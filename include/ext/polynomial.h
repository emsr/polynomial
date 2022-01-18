
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
 * @file polynomial.h Class declaration for a dense univariate polynomial.
 *
 * This file is a GNU extension to the Standard C++ Library.
 *
 * This file contains the declaration of a dense-polynomial class.
 */

/**
 * @def  POLYNOMIAL_H
 *
 * @brief  A guard for the polynomial class header.
 */
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H 1

#include <initializer_list>
#include <vector>
#include <iosfwd>
#include <limits>
#include <array>
#include <utility> // For exchange.
#include <type_traits>
#include <complex>

#include <ext/notsospecfun.h>

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
namespace emsr
{

  template<typename, typename = std::void_t<>>
    struct has_imag_t
    : std::false_type
    { };

  template<typename Tp>
    struct has_imag_t<Tp, std::void_t<decltype(std::declval<Tp&>().imag())>>
    : std::true_type
    { };

  template<typename Tp>
    constexpr auto has_imag_v = has_imag_t<Tp>::value;


  template<typename, typename = std::void_t<>>
    struct has_value_type_t
    : std::false_type
    { };

  template<typename Tp>
    struct has_value_type_t<Tp, std::void_t<typename Tp::value_type>>
    : std::true_type
    { };

  template<typename Tp>
    constexpr auto has_value_type_v = has_value_type_t<Tp>::value;


  template<typename Tp>
    class Polynomial;

  template<typename Tp>
    struct real_type
    { using type = Tp; };

  template<typename Tp>
    struct real_type<std::complex<Tp>>
    { using type = Tp; };

  template<typename Tp>
    struct real_type<Polynomial<Tp>>;

  template<typename Tp>
    using real_type_t = typename real_type<Tp>::type;


  /**
   * @brief A dense polynomial class with a contiguous array of coefficients.
   * The coefficients are lowest-order first:
   * @f[
   *    P(x) = a_0 + a_1 x + ... + a_n x^n
   * @f]
   */
  template<typename Tp>
    class Polynomial
    {
    public:
      /**
       * Typedefs.
       */
      using value_type = typename std::vector<Tp>::value_type;
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
      using real_type = real_type_t<Tp>;

      /**
       * Create a zero degree polynomial with coefficient value zero.
       */
      Polynomial()
      : m_coeff(1)
      { }

      /**
       * Copy ctor.
       */
      Polynomial(const Polynomial&) = default;

      /**
       * Move ctor.
       */
      Polynomial(Polynomial&&) noexcept = default;

      template<typename Up>
	Polynomial(const Polynomial<Up>& poly)
	: m_coeff{}
	{
          for (const auto c : poly)
	    this->m_coeff.push_back(static_cast<value_type>(c));
          this->m_set_scale();
	}

      /**
       * Create a monomial.
       */
      explicit
      Polynomial(value_type a, size_type degree = 0)
      : m_coeff(degree + 1)
      { this->m_coeff[degree] = a; }

      /**
       * Create a polynomial from an initializer list of coefficients.
       */
      Polynomial(std::initializer_list<value_type> ila)
      : m_coeff(ila)
      { this->m_set_scale(); }

      /**
       * Create a polynomial from an input iterator range of coefficients.
       */
      template<typename InIter,
	       typename = std::_RequireInputIter<InIter>>
	Polynomial(const InIter& abegin, const InIter& aend)
	: m_coeff(abegin, aend)
	{ this->m_set_scale(); }

      /**
       * Use Lagrange interpolation to construct a polynomial passing through
       * the data points.  The degree will be one less than the number of points.
       */
      template<typename InIter,
	       typename = std::_RequireInputIter<InIter>>
	Polynomial(const InIter& xbegin, const InIter& xend,
		    const InIter& ybegin)
	: m_coeff()
	{
	  std::vector<Polynomial<value_type>> numer;
	  std::vector<Polynomial<value_type>> denom;
	  for (auto xi = xbegin; xi != xend; ++xi)
	    {
	      for (auto xj = xi + 1; xj != xend; ++xj)
		denom.push_back(value_type(*xj) - value_type(*xi));
	      numer.push_back({-value_type(*xi), value_type{1}});
	    }
          this->m_set_scale();
	}

      /**
       * Create a polynomial from a generator and a maximum degree.
       */
      template<typename Gen>
	Polynomial(Gen gen, size_type degree)
	: m_coeff()
	{
	  this->m_coeff.reserve(degree + 1);
	  for (size_type k = 0; k <= degree; ++k)
	    this->m_coeff.push_back(gen(k));
          this->m_set_scale();
	}

      /**
       * Swap the polynomial with another polynomial.
       */
      void
      swap(Polynomial& poly) noexcept
      { this->m_coeff.swap(poly.m_coeff); }

      /**
       * Evaluate the polynomial at the input point.
       */
      value_type
      operator()(value_type x) const
      {
	if (this->degree() > 0)
	  {
	    if (std::abs(x) <= real_type{1})
	      {
		value_type poly(this->coefficient(this->degree()));
		for (int i = this->degree() - 1; i >= 0; --i)
		  poly = poly * x + this->coefficient(i);
		return poly;
	      }
	    else
	      {
		const auto rx = real_type{1} / x;
		value_type poly(this->coefficient(0));
		for (std::size_t i = 1ull; i <= this->degree(); ++i)
		  poly = poly * rx + this->coefficient(i);
		for (std::size_t i = 1ull; i <= this->degree(); ++i)
		  poly *= x;
		return poly;
	      }
	  }
	else
	  return this->coefficient(0);
      }

      /**
       * Evaluate the polynomial at the input point.
       */
      template<typename Up>
	auto
	operator()(Up x) const
	-> decltype(value_type{} * Up{})
	{
	  if (this->degree() > 0)
	    {
	      if (std::abs(x) <= real_type{1})
		{
		  auto poly(Up{1} * this->coefficient(this->degree()));
		  for (int i = this->degree() - 1; i >= 0; --i)
		    poly = poly * x + this->coefficient(i);
		  return poly;
		}
	      else
		{
		  const auto rx = real_type{1} / x;
		  auto poly(Up{1} * this->coefficient(0));
		  for (std::size_t i = 1ull; i <= this->degree(); ++i)
		    poly = poly * rx + this->coefficient(i);
		  for (std::size_t i = 1ull; i <= this->degree(); ++i)
		    poly *= x;
		  return poly;
		}
	    }
	  else
	    return Up{1} * this->coefficient(0);
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
      template<typename Up>
	auto
	operator()(const std::complex<Up>& z) const
	-> std::enable_if_t<!has_imag_v<Tp>,
			    std::complex<std::decay_t<
		decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>;

      /**
       * Evaluate the polynomial at a range of input points.
       * The output is written to the output iterator which
       * must be large enough to contain the results.
       * The next available output iterator is returned.
       */
      template<typename InIter, typename OutIter,
	       typename = std::_RequireInputIter<InIter>>
	OutIter
	operator()(const InIter& xbegin, const InIter& xend,
        	   OutIter& pbegin) const
	{
	  for (; xbegin != xend; ++xbegin)
	    pbegin++ = (*this)(xbegin++);
	  return pbegin;
	}

      template<size_type N>
	void
	eval(value_type x, std::array<value_type, N>& arr);

      /**
       * Evaluate the polynomial and its derivatives at the point x.
       * The values are placed in the output range starting with the
       * polynomial value and continuing through higher derivatives.
       */
      template<typename OutIter>
	void
	eval(value_type x, OutIter b, OutIter e);

      /**
       * Evaluate the even part of the polynomial at the input point.
       */
      value_type
      eval_even(value_type x) const;

      /**
       * Evaluate the odd part of the polynomial at the input point.
       */
      value_type
      eval_odd(value_type x) const;

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
      template<typename Up>
	auto
	eval_even(const std::complex<Up>& z) const
	-> std::enable_if_t<!has_imag_v<Tp>,
			    std::complex<std::decay_t<
		decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>;

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
      template<typename Up>
	auto
	eval_odd(const std::complex<Up>& z) const
	-> std::enable_if_t<!has_imag_v<Tp>,
			std::complex<std::decay_t<
		decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>;

      /**
       * Return the derivative polynomial.
       */
      Polynomial
      derivative() const
      {
	Polynomial res(value_type{},
			  this->degree() > 0UL ? this->degree() - 1 : 0UL);
	for (size_type n = this->degree(), i = 1; i <= n; ++i)
	  res.m_coeff[i - 1] = i * this->m_coeff[i];
	return res;
      }

      /**
       * Return the derivative of the polynomial at the given point.
       */
      template<typename Up>
        decltype(Up{} * value_type{})
        derivative(Up c) const
        {
	  using res_t = decltype(Up{} * value_type{});
	  const int n = this->degree();
	  res_t res = real_type(n) * this->m_coeff[n];
          for (int i = n - 1; i > 0; --i)
	    res = c * res + real_type(i) * this->m_coeff[i];
	  return res;
        }

      /**
       * Return the integral polynomial with given integration constant.
       */
      Polynomial
      integral(value_type c = value_type{}) const
      {
	Polynomial res(value_type{}, this->degree() + 1);
	res.m_coeff[0] = c;
	for (size_type n = this->degree(), i = 0; i <= n; ++i)
	  res.m_coeff[i + 1] = this->m_coeff[i] / value_type(i + 1);
	return res;
      }

      /**
       * Return the integral of the polynomial with given integration limits.
       */
      template<typename Up>
        decltype(Up{} * value_type{})
	integral(Up a, Up b) const
        {
	  using res_t = decltype(Up{} * value_type{});
	  const int n = this->degree();
	  const auto coeff = this->m_coeff[n] / real_type(n + 1);
	  res_t resa = coeff * a;
	  res_t resb = coeff * b;
	  for (int i = n - 1; i >= 0; --i)
	    {
	      const auto coeff = this->m_coeff[i] / real_type(i + 1);
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
      Polynomial
      operator+() const noexcept
      { return *this; }

      /**
       * Unary minus.
       */
      Polynomial
      operator-() const
      { return Polynomial(*this) *= value_type(-1); }

      /**
       * Assign from a scalar.
       * The result is a zero degree polynomial equal to the scalar.
       */
      Polynomial&
      operator=(const value_type& x)
      {
	this->m_coeff = {x};
	return *this;
      }

      /**
       * Copy assignment.
       */
      Polynomial&
      operator=(const Polynomial&) = default;

      template<typename Up>
	Polynomial&
	operator=(const Polynomial<Up>& poly)
	{
	  if (&poly != this)
	    {
	      this->m_coeff.clear();
	      for (const auto c : poly)
		this->m_coeff.push_back(static_cast<value_type>(c));
	      return *this;
	    }
	}

      /**
       * Assign from an initialiser list.
       */
      Polynomial&
      operator=(std::initializer_list<value_type> ila)
      {
	this->m_coeff = ila;
	return *this;
      }

      /**
       * Add a scalar to the polynomial.
       */
      template<typename Up>
	Polynomial&
	operator+=(const Up& x)
	{
	  this->m_coeff[0] += static_cast<value_type>(x);
	  return *this;
	}

      /**
       * Subtract a scalar from the polynomial.
       */
      template<typename Up>
	Polynomial&
	operator-=(const Up& x)
	{
	  this->m_coeff[0] -= static_cast<value_type>(x);
	  return *this;
	}

      /**
       * Multiply the polynomial by a scalar.
       */
      template<typename Up>
	Polynomial&
	operator*=(const Up& c)
	{
	  for (size_type i = 0; i < this->m_coeff.size(); ++i)
	    this->m_coeff[i] *= static_cast<value_type>(c);
	  return *this;
	}

      /**
       * Divide the polynomial by a scalar.
       */
      template<typename Up>
	Polynomial&
	operator/=(const Up& c)
	{
	  for (size_type i = 0; i < this->m_coeff.size(); ++i)
	    this->m_coeff[i] /= static_cast<value_type>(c);
	  return *this;
	}

      /**
       * Take the modulus of the polynomial relative to a scalar.
       * The result is always a zero polunomial.
       */
      template<typename Up>
	Polynomial&
	operator%=(const Up&)
	{
	  this->degree(0UL); // Resize.
	  this->m_coeff[0] = value_type{};
	  return *this;
	}

      /**
       * Add another polynomial to the polynomial.
       */
      template<typename Up>
	Polynomial&
	operator+=(const Polynomial<Up>& poly)
	{
	  this->degree(std::max(this->degree(), poly.degree()));
	  for (size_type n = poly.degree(), i = 0; i <= n; ++i)
	    this->m_coeff[i] += static_cast<value_type>(poly.m_coeff[i]);
	  return *this;
	}

      /**
       * Subtract another polynomial from the polynomial.
       */
      template<typename Up>
	Polynomial&
	operator-=(const Polynomial<Up>& poly)
	{
	  // Resize if necessary.
	  this->degree(std::max(this->degree(), poly.degree()));
	  for (size_type n = poly.degree(), i = 0; i <= n; ++i)
	    this->m_coeff[i] -= static_cast<value_type>(poly.m_coeff[i]);
	  return *this;
	}

      /**
       * Multiply the polynomial by another polynomial.
       */
      template<typename Up>
	Polynomial&
	operator*=(const Polynomial<Up>& poly);

      /**
       * Divide the polynomial by another polynomial.
       */
      template<typename Up>
	Polynomial&
	operator/=(const Polynomial<Up>& poly)
	{
	  Polynomial<value_type >quo, rem;
	  divmod(*this, poly, quo, rem);
	  *this = quo;
	  return *this;
	}

      /**
       * Take the modulus of (modulate?) the polynomial relative to another polynomial.
       */
      template<typename Up>
	Polynomial&
	operator%=(const Polynomial<Up>& poly)
	{
	  Polynomial<value_type >quo, rem;
	  divmod(*this, poly, quo, rem);
	  *this = rem;
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
	    this->m_coeff[n - i] += shift * this->m_coeff[n - i + 1];
      }

      /**
       * Return the degree or the power of the largest coefficient.
       */
      size_type
      degree() const noexcept
      { return this->m_coeff.size() - 1; }

      /**
       * Set the degree or the power of the largest coefficient.
       */
      void
      degree(size_type degree)
      { this->m_coeff.resize(degree + 1UL); }

      /**
       * Return the size of the coefficient sequence.
       */
      size_type
      size() const noexcept
      { return this->m_coeff.size(); }

      /**
       * Return the @c ith coefficient with range checking.
       */
      value_type
      coefficient(size_type i) const
      { return this->m_coeff.at(i); }

      /**
       * Set coefficient @c i to @c val with range checking.
       */
      void
      coefficient(size_type i, value_type val)
      { this->m_coeff.at(i) = val; }

      /**
       * Return coefficient @c i.
       */
      value_type
      operator[](size_type i) const noexcept
      { return this->m_coeff[i]; }

      /**
       * Return coefficient @c i as an assignable quantity.
       */
      reference
      operator[](size_type i) noexcept
      { return this->m_coeff[i]; }

      /**
       * Return a const vector of coefficients.
       */
      const std::vector<value_type>
      coefficients() const noexcept
      { return this->m_coeff; }

      /**
       * Return a vector of coefficients.
       */
      std::vector<value_type>
      coefficients() noexcept
      { return this->m_coeff; }

      /**
       * Return a @c const pointer to the coefficient sequence.
       */
      const value_type*
      data() const noexcept
      { return this->m_coeff.data(); }

      /**
       * Return a @c pointer to the coefficient sequence.
       */
      value_type*
      data() noexcept
      { return this->m_coeff.data(); }

      /**
       * Return an iterator to the beginning of the coefficient sequence.
       */
      iterator
      begin() noexcept
      { return this->m_coeff.begin(); }

      /**
       * Return an iterator to one past the end of the coefficient sequence.
       */
      iterator
      end() noexcept
      { return this->m_coeff.end(); }

      /**
       * Return a @c const iterator the beginning
       * of the coefficient sequence.
       */
      const_iterator
      begin() const noexcept
      { return this->m_coeff.begin(); }

      /**
       * Return a @c const iterator to one past the end
       * of the coefficient sequence.
       */
      const_iterator
      end() const noexcept
      { return this->m_coeff.end(); }

      /**
       * Return a @c const iterator the beginning
       * of the coefficient sequence.
       */
      const_iterator
      cbegin() const noexcept
      { return this->m_coeff.cbegin(); }

      /**
       * Return a @c const iterator to one past the end
       * of the coefficient sequence.
       */
      const_iterator
      cend() const noexcept
      { return this->m_coeff.cend(); }

      reverse_iterator
      rbegin() noexcept
      { return this->m_coeff.rbegin(); }

      reverse_iterator
      rend() noexcept
      { return this->m_coeff.rend(); }

      const_reverse_iterator
      rbegin() const noexcept
      { return this->m_coeff.rbegin(); }

      const_reverse_iterator
      rend() const noexcept
      { return this->m_coeff.rend(); }

      const_reverse_iterator
      crbegin() const noexcept
      { return this->m_coeff.crbegin(); }

      const_reverse_iterator
      crend() const noexcept
      { return this->m_coeff.crend(); }

      template<typename CharT, typename Traits, typename Tp1>
	friend std::basic_istream<CharT, Traits>&
	operator>>(std::basic_istream<CharT, Traits>&, Polynomial<Tp1>&);

      template<typename Tp1>
	friend bool
	operator==(const Polynomial<Tp1>& pa,
		   const Polynomial<Tp1>& pb);

      /**
       * Remove zero max-order coefficients.
       */
      Polynomial&
      deflate(real_type max_abs_coef)
      {
	size_type n = this->degree();
	for (size_type i = this->degree(); i > 0; --i)
	  if (std::abs(this->m_coeff[i]) < max_abs_coef)
	    --n;
	  else
	    break;
	this->degree(n);
        return *this;
      }

      /**
       * Divide the polynomial by an input polynomia and remove zero
       * max-order coefficients.
       */
      Polynomial&
      deflate(const Polynomial<value_type>& poly,
	      real_type max_abs_coef)
      {
	Polynomial<value_type> quo, rem;
	divmod(*this, poly, quo, rem);

	// Remainder should be null.
	size_type n = rem.degree();
	for (size_type i = rem.degree(); i > 0; --i)
	  if (std::abs(rem[i]) < max_abs_coef)
	    --n;
	  else
	    break;

	if (n == 0)
	  *this = quo.deflate(max_abs_coef);
	else
	  throw std::runtime_error("deflate: ");

        return *this;
      }

    private:

      /// Return the scale.
      real_type
      m_get_scale() const
      { return this->m_scale; }

      real_type m_scale = real_type{1};

      void m_set_scale();

      std::vector<value_type> m_coeff;
    };

  template<typename Tp>
    struct real_type<Polynomial<Tp>>
    { using type = typename Polynomial<Tp>::real_type; };

  /**
   * Return the scale for a polynomial.
   */
  template<typename Tp>
    real_type_t<Polynomial<Tp>>
    get_scale(const Polynomial<Tp>& poly)
    { return poly.m_get_scale(); }

  /**
   * Return the scale for a number.
   */
  template<typename Tp>
    decltype(std::abs(Tp()))
    get_scale(const Tp& x)
    { return std::abs(x); }

  /**
   * Return the sum of a polynomial with a scalar.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() + Up())>
    operator+(const Polynomial<Tp>& poly, const Up& x)
    { return Polynomial<decltype(Tp() + Up())>(poly) += x; }

  /**
   *
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() + Up())>
    operator+(const Tp& x, const Polynomial<Up>& poly)
    { return Polynomial<decltype(Tp() + Up())>(x) += poly; }

  /**
   * Return the difference of a polynomial with a scalar.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() - Up())>
    operator-(const Polynomial<Tp>& poly, const Up& x)
    { return Polynomial<decltype(Tp() - Up())>(poly) -= x; }

  /**
   *
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() - Up())>
    operator-(const Tp& x, const Polynomial<Up>& poly)
    { return Polynomial<decltype(Tp() - Up())>(x) -= poly; }

  /**
   * Return the product of a polynomial with a scalar.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() * Up())>
    operator*(const Polynomial<Tp>& poly, const Up& x)
    { return Polynomial<decltype(Tp() * Up())>(poly) *= x; }

  /**
   *
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() * Up())>
    operator*(const Tp& x, const Polynomial<Up>& poly)
    { return Polynomial<decltype(Tp() * Up())>(x) *= poly; }

  /**
   * Return the quotient of a polynomial with a scalar.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() / Up())>
    operator/(const Polynomial<Tp>& poly, const Up& x)
    { return Polynomial<decltype(Tp() / Up())>(poly) /= x; }

  /**
   *
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() / Up())>
    operator%(const Polynomial<Tp>& poly, const Up& x)
    { return Polynomial<decltype(Tp() / Up())>(poly) %= x; }

  /**
   * Return the sum of two polynomials.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() + Up())>
    operator+(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
    { return Polynomial<decltype(Tp() + Up())>(pa) += pb; }

  /**
   * Return the difference of two polynomials.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() - Up())>
    operator-(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
    { return Polynomial<decltype(Tp() - Up())>(pa) -= pb; }

  /**
   * Return the product of two polynomials.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() * Up())>
    operator*(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
    { return Polynomial<decltype(Tp() * Up())>(pa) *= pb; }

  /**
   * Return the quotient of two polynomials.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() / Up())>
    operator/(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
    { return Polynomial<decltype(Tp() / Up())>(pa) /= pb; }

  /**
   * Return the modulus or remainder of one polynomial relative to another one.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() / Up())>
    operator%(const Polynomial<Tp>& pa, const Polynomial<Up>& pb)
    { return Polynomial<decltype(Tp() / Up())>(pa) %= pb; }

  /**
   * Return the quotient of a scalar and a polynomials.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() / Up())>
    operator/(const Tp& x, const Polynomial<Up>& poly)
    { return Polynomial<decltype(Tp() / Up())>(x) /= poly; }

  /**
   * Return the modulus or remainder of a scalar divided by a polynomial.
   */
  template<typename Tp, typename Up>
    inline Polynomial<decltype(Tp() / Up())>
    operator%(const Tp& x, const Polynomial<Up>& poly)
    { return Polynomial<decltype(Tp() / Up())>(x) %= poly; }

  /**
   * Divide two polynomials returning the quotient and remainder.
   */
  template<typename Tp>
    void
    divmod(const Polynomial<Tp>& num, const Polynomial<Tp>& den,
           Polynomial<Tp>& quo, Polynomial<Tp>& rem);

  /**
   * Write a polynomial to a stream.
   * The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename Tp>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os,
	       const Polynomial<Tp>& poly);

  /**
   * Read a polynomial from a stream.
   * The input format can be a plain scalar (zero degree polynomial)
   * or a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename Tp>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& is,
	       Polynomial<Tp>& poly);

  /**
   * Return true if two polynomials are equal.
   */
  template<typename Tp>
    inline bool
    operator==(const Polynomial<Tp>& pa, const Polynomial<Tp>& pb)
    { return pa.m_coeff == pb.m_coeff; }

  /**
   * Return false if two polynomials are equal.
   */
  template<typename Tp>
    inline bool
    operator!=(const Polynomial<Tp>& pa, const Polynomial<Tp>& pb)
    { return !(pa == pb); }

  /**
   * See Polynomial::swap().
   */
  template<typename Tp>
    inline void
    swap(Polynomial<Tp>& pa, Polynomial<Tp>& pb)
    noexcept(noexcept(pa.swap(pb)))
    { pa.swap(pb); }

} // namespace emsr

#include <ext/polynomial.tcc>

#endif // POLYNOMIAL_H

