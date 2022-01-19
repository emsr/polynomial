
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
 * @file solution.h
 *
 * This file contains a type representing a solution of a polynomial.
 * The soution could be null, i.e. on-existent.
 * If it exists it could be real of complex.
 */

/**
 * @def  SOLUTION_H
 *
 * @brief  A guard for the polynomial solution class header.
 */
#ifndef SOLUTION_H
#define SOLUTION_H 1

#include <complex>
#include <variant>
#include <iosfwd>

namespace emsr
{

  template<typename Real>
    struct Solution
    : public std::variant<std::monostate, Real, std::complex<Real>>
    {
      using Base = std::variant<std::monostate, Real, std::complex<Real>>;

      Solution() noexcept
      : Base()
      {}

      Solution(Real x) noexcept
      : Base(x)
      {}

      Solution(const std::complex<Real>& x) noexcept
      : Base(x)
      {}
    };

  template<typename Real>
    constexpr bool
    is_valid(const Solution<Real>& x)
    { return x.index() != 0; }

  template<typename Real>
    constexpr Real
    real(const Solution<Real>& x)
    {
      if (x.index() == 0)
	return std::numeric_limits<Real>::quiet_NaN();
      else if (x.index() == 1)
	return std::get<1>(x);
      else
	return std::real(std::get<2>(x));
    }

  template<typename Real>
    constexpr Real
    imag(const Solution<Real>& x)
    {
      if (x.index() == 0)
	return std::numeric_limits<Real>::quiet_NaN();
      else if (x.index() == 1)
	return Real{0};
      else
	return std::imag(std::get<2>(x));
    }

  template<typename Real>
    constexpr Real
    abs(const Solution<Real>& x)
    {
      if (x.index() == 0)
	return std::numeric_limits<Real>::quiet_NaN();
      else if (x.index() == 1)
	return std::abs(std::get<1>(x));
      else
	return std::abs(std::get<2>(x));
    }

  /**
   * Return the solution as a complex number or NaN.
   */
  template<typename Real>
    constexpr Solution<Real>
    to_complex(const Solution<Real>& x)
    {
      if (x.index() == 0)
	return Solution<Real>();
      else if (x.index() == 1)
	return Solution<Real>(std::complex<Real>(std::get<1>(x)));
      else
	return x;
    }

  /**
   * Unary +/-
   */
  template<typename Real>
    constexpr Solution<Real>
    operator+(const Solution<Real>& x)
    { return x; }

  template<typename Real>
    constexpr Solution<Real>
    operator-(const Solution<Real>& x)
    {
      if (x.index() == 0)
	return x;
      else if (x.index() == 1)
	return Solution<Real>(-std::get<1>(x));
      else
	return Solution<Real>(-std::get<2>(x));
    }

  /**
   * Addition operators...
   */
  template<typename Real>
    constexpr Solution<Real>
    operator+(const Solution<Real>& x, const Solution<Real>& y)
    {
      if (x.index() == 0)
	return x;
      if (y.index() == 0)
	return y;
      else if (x.index() == 1)
	{
	  if (y.index() == 1)
	    return Solution<Real>(std::get<1>(x) + std::get<1>(y));
	  else
	    return Solution<Real>(std::get<1>(x) + std::get<2>(y));
	}
      else
	{
	  if (y.index() == 1)
	    return Solution<Real>(std::get<2>(x) + std::get<1>(y));
	  else
	    return Solution<Real>(std::get<2>(x) + std::get<2>(y));
	}
    }

  template<typename Real>
    constexpr Solution<Real>
    operator+(const Solution<Real>& x, Real y)
    { return operator+(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr Solution<Real>
    operator+(Real x, const Solution<Real>& y)
    { return operator+(Solution<Real>(x), y); }

  template<typename Real>
    constexpr Solution<Real>
    operator+(const Solution<Real>& x, const std::complex<Real>& y)
    { return operator+(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr Solution<Real>
    operator+(const std::complex<Real>& x, const Solution<Real>& y)
    { return operator+(Solution<Real>(x), y); }

  /**
   * Subtraction operators...
   */
  template<typename Real>
    constexpr Solution<Real>
    operator-(const Solution<Real>& x, const Solution<Real>& y)
    {
      if (x.index() == 0)
	return x;
      if (y.index() == 0)
	return y;
      else if (x.index() == 1)
	{
	  if (y.index() == 1)
	    return Solution<Real>(std::get<1>(x) - std::get<1>(y));
	  else
	    return Solution<Real>(std::get<1>(x) - std::get<2>(y));
	}
      else
	{
	  if (y.index() == 1)
	    return Solution<Real>(std::get<2>(x) - std::get<1>(y));
	  else
	    return Solution<Real>(std::get<2>(x) - std::get<2>(y));
	}
    }

  template<typename Real>
    constexpr Solution<Real>
    operator-(const Solution<Real>& x, Real y)
    { return operator-(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr Solution<Real>
    operator-(Real x, const Solution<Real>& y)
    { return operator-(Solution<Real>(x), y); }

  template<typename Real>
    constexpr Solution<Real>
    operator-(const Solution<Real>& x, const std::complex<Real>& y)
    { return operator-(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr Solution<Real>
    operator-(const std::complex<Real>& x, const Solution<Real>& y)
    { return operator-(Solution<Real>(x), y); }

  /**
   * Multiplication operators...
   */
  template<typename Real>
    constexpr Solution<Real>
    operator*(const Solution<Real>& x, const Solution<Real>& y)
    {
      if (x.index() == 0)
	return x;
      if (y.index() == 0)
	return y;
      else if (x.index() == 1)
	{
	  if (y.index() == 1)
	    return Solution<Real>(std::get<1>(x) * std::get<1>(y));
	  else
	    return Solution<Real>(std::get<1>(x) * std::get<2>(y));
	}
      else
	{
	  if (y.index() == 1)
	    return Solution<Real>(std::get<2>(x) * std::get<1>(y));
	  else
	    return Solution<Real>(std::get<2>(x) * std::get<2>(y));
	}
    }

  template<typename Real>
    constexpr Solution<Real>
    operator*(const Solution<Real>& x, Real y)
    { return operator*(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr Solution<Real>
    operator*(Real x, const Solution<Real>& y)
    { return operator*(Solution<Real>(x), y); }

  template<typename Real>
    constexpr Solution<Real>
    operator*(const Solution<Real>& x, const std::complex<Real>& y)
    { return operator*(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr Solution<Real>
    operator*(const std::complex<Real>& x, const Solution<Real>& y)
    { return operator*(Solution<Real>(x), y); }

  /**
   * division operators...
   */
  template<typename Real>
    constexpr Solution<Real>
    operator/(const Solution<Real>& x, const Solution<Real>& y)
    {
      if (x.index() == 0)
	return x;
      if (y.index() == 0)
	return y;
      else if (x.index() == 1)
	{
	  if (y.index() == 1)
	    return Solution<Real>(std::get<1>(x) / std::get<1>(y));
	  else
	    return Solution<Real>(std::get<1>(x) / std::get<2>(y));
	}
      else
	{
	  if (y.index() == 1)
	    return Solution<Real>(std::get<2>(x) / std::get<1>(y));
	  else
	    return Solution<Real>(std::get<2>(x) / std::get<2>(y));
	}
    }

  template<typename Real>
    constexpr Solution<Real>
    operator/(const Solution<Real>& x, Real y)
    { return operator/(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr Solution<Real>
    operator/(Real x, const Solution<Real>& y)
    { return operator/(Solution<Real>(x), y); }

  template<typename Real>
    constexpr Solution<Real>
    operator/(const Solution<Real>& x, const std::complex<Real>& y)
    { return operator/(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr Solution<Real>
    operator/(const std::complex<Real>& x, const Solution<Real>& y)
    { return operator/(Solution<Real>(x), y); }

  /**
   * Test for equality and inequality.
   */
  template<typename Real>
    constexpr bool
    operator==(const Solution<Real>& x, const Solution<Real>& y)
    {
      if (x.index() == 0 || y.index() == 0)
	return false;
      else if (x.index() == y.index())
	{
	  if (x.index() == 1)
	    return std::get<1>(x) == std::get<1>(y);
	  else
	    return std::get<2>(x) == std::get<2>(y);
	}
      else
	return false;
    }

  template<typename Real>
    bool
    operator==(const Solution<Real>& x, Real y)
    { return x == Solution<Real>(y); }

  template<typename Real>
    constexpr bool
    operator==(Real x, const Solution<Real>& y)
    { return Solution<Real>(x) == y; }

  template<typename Real>
    constexpr bool
    operator==(const Solution<Real>& x, const std::complex<Real>& y)
    { return x == Solution<Real>(y); }

  template<typename Real>
    constexpr bool
    operator==(const std::complex<Real>& x, const Solution<Real>& y)
    { return Solution<Real>(x) == y; }

  template<typename Real>
    constexpr bool
    operator!=(const Solution<Real>& x, const Solution<Real>& y)
    { return !(x == y); }

  template<typename Real>
    bool
    operator!=(const Solution<Real>& x, Real y)
    { return !(x == y); }

  template<typename Real>
    constexpr bool
    operator!=(Real x, const Solution<Real>& y)
    { return !(x == y); }

  template<typename Real>
    constexpr bool
    operator!=(const Solution<Real>& x, const std::complex<Real>& y)
    { return !(x == y); }

  template<typename Real>
    constexpr bool
    operator!=(const std::complex<Real>& x, const Solution<Real>& y)
    { return !(x == y); }

  /**
   * Lexicographic order of solutions as complex numbers.
   * Null solutions compare as less than except to another null solution.
   *
   * A tribool might be a good thing for this when either
   * of the solutions is null.
   */
  template<typename Real>
    constexpr bool
    operator<(const Solution<Real>& x, const Solution<Real>& y)
    {
      if (x.index() == 0 && y.index() == 0)
	return false;
      else if (x.index() == 0)
	return true;
      else if (y.index() == 0)
	return false;
      else
	{
	  const auto rex = emsr::real(x);
	  const auto rey = emsr::real(y);
	  if (rex < rey)
	    return true;
	  else if (rex == rey)
	    return emsr::imag(x) < emsr::imag(y);
	  else
	    return false;
	}
    }

  template<typename Real>
    constexpr bool
    operator<(const Solution<Real>& x, Real y)
    { return operator<(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr bool
    operator<(Real x, const Solution<Real>& y)
    { return operator<(Solution<Real>(x), y); }

  template<typename Real>
    constexpr bool
    operator<(const Solution<Real>& x, const std::complex<Real>& y)
    { return operator<(x, Solution<Real>(y)); }

  template<typename Real>
    constexpr bool
    operator<(const std::complex<Real>& x, const Solution<Real>& y)
    { return operator<(Solution<Real>(x), y); }

  /**
   * Output a solution to a stream.
   */
  template<typename Char, typename Real>
    std::basic_ostream<Char>&
    operator<<(std::basic_ostream<Char>& out,
	       const emsr::Solution<Real>& sln)
    {
      const auto idx = sln.index();
      if (idx == 0)
	out << "null";
      else if (idx == 1)
	out << std::get<1>(sln);
      else if (idx == 2)
	out << std::get<2>(sln);
      return out;
    }

} // namespace emsr

#endif // SOLUTION_H
