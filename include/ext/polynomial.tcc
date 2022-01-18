
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
 * @file polynomial.tcc Out-of-line definitions of members for
 * a dense univariate polynomial.
 *
 * This file is a GNU extension to the Standard C++ Library.
 * This file contains the out-of-line implementations of the polynomial class.
 *
 * @see polynomial.h
 */

/**
 * @def  POLYNOMIAL_TCC
 *
 * @brief  A guard for the polynomial class implementation header.
 */
#ifndef POLYNOMIAL_TCC
#define POLYNOMIAL_TCC 1

#include <ios>
#include <complex>
#include <utility> // For exchange.

namespace emsr
{

// We also need an integer coef version :-\ ?
  template<typename Tp>
    void
    Polynomial<Tp>::m_set_scale()
    { }
  /**
   * 
OTOH, What does it mean to get roots of integer polynomials - the roots will
often be real or even complex.

What does it even mean to get roots of polynomial polynomials.
I think you'd want to pick a parameter for the inner coefficient polynomials.

Wait! This has nothing to do with the coefficient type.  It's the evaluation type
that determines this. It may be sane to make real_type be double if it comes in
as integral.

OTOOH, scaling down an integer polynomial will either screw up the polynomial
or it will make a polynomial of a different type.

This scaling thing can only apply to real or complex polynomials.

  template<typename Tp>
    std::enable_if_t<std::is_integral_v<Tp>>
    Polynomial<Tp>::m_set_scale()
    { }
   */

  /**
   * 
  template<typename Tp>
    //std::enable_if_t<std::is_floating_point_v<Tp>>
    void
    Polynomial<Tp>::m_set_scale()
    {
      constexpr real_type s_eps = std::numeric_limits<real_type>::epsilon();
      constexpr real_type
	s_base = real_type{std::numeric_limits<real_type>::radix};
      constexpr real_type s_tiny = s_eps * s_eps * s_eps; // ~ 1.0e-50;
      constexpr real_type s_huge = std::numeric_limits<real_type>::max();
      constexpr real_type s_low = s_tiny / s_eps;

      // Find largest and smallest moduli of coefficients.
      auto a_max = real_type{0};
      auto a_min = s_huge;
      for (int i = 0; i <= this->degree(); ++i)
	{
	  const auto x = get_scale(this->m_coeff[i]);
	  if (x > a_max)
	    a_max = x;
	  if (x != real_type{0} && x < a_min)
	    a_min = x;
	}
      // Scale if there are large or very tiny coefficients.
      // Computes a scale factor to multiply the coefficients
      // of the polynomial. The scaling is done to avoid overflow
      // and to avoid undetected underflow interfering
      // with the convergence criterion.
      // The factor is a power of the base.
      auto scale = s_low / a_min;
      bool rescale = true;
      if (scale > real_type{1} && s_huge / scale < a_max)
	rescale = false;
      if (scale <= real_type{1} && a_max < real_type{10})
	rescale = false;

      this->m_scale = real_type{1};
      if (rescale)
	{
	  // Scale polynomial.
	  if (scale == real_type{0})
	    scale = s_tiny;
	  const auto lg = std::ilogb(scale);
	  this->m_scale = std::pow(s_base, lg);
	  //if (this->m_scale != real_type{1})
	  //  for (int i = 0; i <= this->degree(); ++i)
	  //    this->m_coeff[i] *= this->m_scale;
	}
    }
   */

  /**
   * Scaling specialization for polynomial coefficient polynomials.
  template<typename Tp>
    void
    Polynomial<Polynomial<Tp>>::m_set_scale()
    {
      for (int i = 0; i <= this->degree(); ++i)
	this->m_coeff[i].m_set_scale();
    }
   */

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
  template<typename Tp>
    template<typename Up>
      auto
      Polynomial<Tp>::operator()(const std::complex<Up>& z) const
      -> std::enable_if_t<!has_imag_v<Tp>,
			  std::complex<std::decay_t<
		decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>
      {
	const auto r = Tp{2} * std::real(z);
	const auto s = std::norm(z);
	size_type n = this->degree();
	auto aa = this->coefficient(n);
	if (n > 0)
	  {
	    auto bb = this->coefficient(n - 1);
	    for (size_type j = 2; j <= n; ++j)
	      bb = this->coefficient(n - j)
		   - s * std::exchange(aa, bb + r * aa);
	    return aa * z + bb;
	  }
	else
	  return aa;
      };

  //  Could/should this be done by output iterator range?
  template<typename Tp>
    template<typename Polynomial<Tp>::size_type N>
      void
      Polynomial<Tp>::eval(typename Polynomial<Tp>::value_type x,
			     std::array<Polynomial<Tp>::value_type, N>& arr)
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
		  arr[j] = std::fma(arr[j], x, arr[j - 1]);
		arr[0] = std::fma(arr[0], x, this->coefficient(i));
	      }
	    //  Now put in the factorials.
	    value_type fact = value_type(1);
	    for (size_type n = arr.size(), i = 2; i < n; ++i)
	      {
		fact *= value_type(i);
		arr[i] *= fact;
	      }
	  }
      }

  /**
   * Evaluate the polynomial and its derivatives at the point x.
   * The values are placed in the output range starting with the
   * polynomial value and continuing through higher derivatives.
   */
  template<typename Tp>
    template<typename OutIter>
      void
      Polynomial<Tp>::eval(typename Polynomial<Tp>::value_type x,
			     OutIter b, OutIter e)
      {
	if(b != e)
	  {
	    std::fill(b, e, value_type{});
	    const size_type sz = m_coeff.size();
	    *b = m_coeff[sz - 1];
            for (int i = sz - 2; i >= 0; --i)
	      {
		for (auto it = std::reverse_iterator<OutIter>(e);
		     it != std::reverse_iterator<OutIter>(b) - 1; ++it)
		  *it = std::fma(*it, x, *(it + 1));
		*b = std::fma(*b, x, m_coeff[i]);
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
   * Evaluate the even part of the polynomial at the input point.
   */
  template<typename Tp>
    typename Polynomial<Tp>::value_type
    Polynomial<Tp>::eval_even(typename Polynomial<Tp>::value_type x) const
    {
      if (this->degree() > 0)
	{
	  const auto odd = this->degree() % 2;
	  const auto xx = x * x;
	  auto poly(this->coefficient(this->degree() - odd));
	  for (int i = this->degree() - odd - 2; i >= 0; i -= 2)
	    poly = std::fma(xx, poly, this->coefficient(i));
	  return poly;
	}
      else
	return value_type{};
    }

  /**
   * Evaluate the odd part of the polynomial at the input point.
   */
  template<typename Tp>
    typename Polynomial<Tp>::value_type
    Polynomial<Tp>::eval_odd(typename Polynomial<Tp>::value_type x) const
    {
      if (this->degree() > 0)
	{
	  const auto even = (this->degree() % 2 == 0 ? 1 : 0);
	  const auto xx = x * x;
	  auto poly(this->coefficient(this->degree() - even));
	  for (int i = this->degree() - even - 2; i >= 0; i -= 2)
	    poly = std::fma(xx, poly, this->coefficient(i));
	  return x * poly;
	}
      else
	return value_type{};
    }

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
  template<typename Tp>
    template<typename Up>
      auto
      Polynomial<Tp>::eval_even(const std::complex<Up>& z) const
      -> std::enable_if_t<!has_imag_v<Tp>,
			  std::complex<std::decay_t<
		decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>
      {
	using real_t = std::decay_t<decltype(value_type{} * Up{})>;
	using cmplx_t = std::complex<real_t>;
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
	      bb = std::fma(-s, std::exchange(aa, bb + r * aa),
			      this->coefficient(n - j));
	    return std::fma(cmplx_t(aa), cmplx_t(zz), cmplx_t(bb));
	  }
	else
	  return cmplx_t{};
      };

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
  template<typename Tp>
    template<typename Up>
      auto
      Polynomial<Tp>::eval_odd(const std::complex<Up>& z) const
      -> std::enable_if_t<!has_imag_v<Tp>,
			  std::complex<std::decay_t<
		decltype(typename Polynomial<Tp>::value_type{} * Up{})>>>
      {
	using real_t = std::decay_t<decltype(value_type{} * Up{})>;
	using cmplx_t = std::complex<real_t>;
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
	      bb = std::fma(-s, std::exchange(aa, bb + r * aa),
			      this->coefficient(n - j));
	    return z
		 * std::fma(cmplx_t(aa), cmplx_t(zz), cmplx_t(bb));
	  }
	else
	  return cmplx_t{};
      };

    /**
     * Multiply the polynomial by another polynomial.
     */
  template<typename Tp>
    template<typename Up>
      Polynomial<Tp>&
      Polynomial<Tp>::operator*=(const Polynomial<Up>& poly)
      {
	//  Test for zero size polys and do special processing?
	const size_type m = this->degree();
	const size_type n = poly.degree();
	std::vector<value_type> new_coeff(m + n + 1);
	for (size_type i = 0; i <= m; ++i)
	  for (size_type j = 0; j <= n; ++j)
	    new_coeff[i + j] += this->m_coeff[i]
				* static_cast<value_type>(poly.m_coeff[j]);
	this->m_coeff = new_coeff;
	return *this;
      }

  /**
   * Divide two polynomials returning the quotient and remainder.
   */
  template<typename Tp>
    void
    divmod(const Polynomial<Tp>& num, const Polynomial<Tp>& den,
           Polynomial<Tp>& quo, Polynomial<Tp>& rem)
    {
      rem = num;
      quo = Polynomial<Tp>(Tp{}, num.degree());
      const std::size_t d_num = num.degree();
      const std::size_t d_den = den.degree();
      if (d_den <= d_num)
	{
	  for (int k = d_num - d_den; k >= 0; --k)
	    {
	      quo.coefficient(k, rem.coefficient(d_den + k)
				   / den.coefficient(d_den));
	      for (int j = d_den + k - 1; j >= k; --j)
		rem.coefficient(j, rem.coefficient(j)
				     - quo.coefficient(k)
				     * den.coefficient(j - k));
	    }
	  quo.degree(d_num - d_den);
	  rem.degree(d_den > 0 ? d_den - 1 : 0);
	}
      else
	quo.degree(0);
    }

  /**
   * Write a polynomial to a stream.
   * The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename Tp>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os,
	       const Polynomial<Tp>& poly)
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
   * Read a polynomial from a stream.
   * The input format can be a plain scalar (zero degree polynomial)
   * or a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename Tp>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& is,
	       Polynomial<Tp>& poly)
    {
      Tp x;
      CharT ch;
      is >> ch;
      if (ch == '(')
	{
	  do
	    {
	      is >> x >> ch;
	      poly.m_coeff.push_back(x);
	    }
	  while (ch == ',');
	  if (ch != ')')
	    is.setstate(std::ios_base::failbit);
	}
      else
	{
	  is.putback(ch);
	  is >> x;
	  poly = x;
	}
      // No null polynomial.
      if (poly.size() == 0)
	poly.m_coeff.resize(1, Tp{});

      return is;
    }

} // namespace emsr

#endif // POLYNOMIAL_TCC

