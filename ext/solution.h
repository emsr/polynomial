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

#ifndef SOLUTION_H
#define SOLUTION_H 1

/**
 * @file solution.h
 *
 * This file is a GNU extension to the Standard C++ Library.
 *
 * This file contains a type representing a solution of a polynomial.
 * The soution could be null, i.e. on-existent.
 * If it exists it could be real of complex.
 */

#pragma GCC system_header

#if __cplusplus < 201703L
# include <bits/c++0x_warning.h>
#else

#include <complex>
#include <variant>
#include <iosfwd>

namespace __gnu_cxx
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Real>
    using solution_t
	    = std::variant<std::monostate, _Real, std::complex<_Real>>;

  template<typename _Real>
    constexpr bool
    is_valid(const solution_t<_Real>& __x)
    { return __x.index() != 0; }
/*
  template<typename _Real>
    constexpr solution_t<_Real>
    real(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return solution_t<_Real>();
      else if (__x.index() == 1)
	return solution_t<_Real>(std::get<1>(__x));
      else
	return solution_t<_Real>(std::real(std::get<2>(__x)));
    }

  template<typename _Real>
    constexpr solution_t<_Real>
    imag(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return solution_t<_Real>();
      else if (__x.index() == 1)
	return solution_t<_Real>(_Real{0});
      else
	return solution_t<_Real>(std::imag(std::get<2>(__x)));
    }
*/
  template<typename _Real>
    constexpr _Real
    real(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return std::numeric_limits<_Real>::quiet_NaN();
      else if (__x.index() == 1)
	return std::get<1>(__x);
      else
	return std::real(std::get<2>(__x));
    }

  template<typename _Real>
    constexpr _Real
    imag(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return std::numeric_limits<_Real>::quiet_NaN();
      else if (__x.index() == 1)
	return _Real{0};
      else
	return std::imag(std::get<2>(__x));
    }

  template<typename _Real>
    constexpr _Real
    abs(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return std::numeric_limits<_Real>::quiet_NaN();
      else if (__x.index() == 1)
	return std::abs(std::get<1>(__x));
      else
	return std::abs(std::get<2>(__x));
    }

  /**
   * Return the solution as a complex number or NaN.
   */
  template<typename _Real>
    constexpr solution_t<_Real>
    to_complex(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return solution_t<_Real>();
      else if (__x.index() == 1)
	return solution_t<_Real>(std::complex<_Real>(std::get<1>(__x)));
      else
	return __x;
    }

  /**
   * Unary +-
   */
  template<typename _Real>
    constexpr solution_t<_Real>
    operator+(const solution_t<_Real>& __x)
    { return __x; }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator-(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return __x;
      else if (__x.index() == 1)
	return solution_t<_Real>(-std::get<1>(__x));
      else
	return solution_t<_Real>(-std::get<2>(__x));
    }

  /**
   * Addition operators...
   */
  template<typename _Real>
    constexpr solution_t<_Real>
    operator+(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      if (__x.index() == 0)
	return __x;
      if (__y.index() == 0)
	return __y;
      else if (__x.index() == 1)
	{
	  if (__y.index() == 1)
	    return solution_t<_Real>(std::get<1>(__x) + std::get<1>(__y));
	  else
	    return solution_t<_Real>(std::get<1>(__x) + std::get<2>(__y));
	}
      else
	{
	  if (__y.index() == 1)
	    return solution_t<_Real>(std::get<2>(__x) + std::get<1>(__y));
	  else
	    return solution_t<_Real>(std::get<2>(__x) + std::get<2>(__y));
	}
    }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator+(const solution_t<_Real>& __x, _Real __y)
    { return operator+(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator+(_Real __x, const solution_t<_Real>& __y)
    { return operator+(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator+(const solution_t<_Real>& __x, std::complex<_Real>& __y)
    { return operator+(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator+(std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator+(solution_t<_Real>(__x), __y); }

  /**
   * Subtraction operators...
   */
  template<typename _Real>
    constexpr solution_t<_Real>
    operator-(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      if (__x.index() == 0)
	return __x;
      if (__y.index() == 0)
	return __y;
      else if (__x.index() == 1)
	{
	  if (__y.index() == 1)
	    return solution_t<_Real>(std::get<1>(__x) - std::get<1>(__y));
	  else
	    return solution_t<_Real>(std::get<1>(__x) - std::get<2>(__y));
	}
      else
	{
	  if (__y.index() == 1)
	    return solution_t<_Real>(std::get<2>(__x) - std::get<1>(__y));
	  else
	    return solution_t<_Real>(std::get<2>(__x) - std::get<2>(__y));
	}
    }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator-(const solution_t<_Real>& __x, _Real __y)
    { return operator-(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator-(_Real __x, const solution_t<_Real>& __y)
    { return operator-(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator-(const solution_t<_Real>& __x, std::complex<_Real>& __y)
    { return operator-(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator-(std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator-(solution_t<_Real>(__x), __y); }

  /**
   * Multiplication operators...
   */
  template<typename _Real>
    constexpr solution_t<_Real>
    operator*(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      if (__x.index() == 0)
	return __x;
      if (__y.index() == 0)
	return __y;
      else if (__x.index() == 1)
	{
	  if (__y.index() == 1)
	    return solution_t<_Real>(std::get<1>(__x) * std::get<1>(__y));
	  else
	    return solution_t<_Real>(std::get<1>(__x) * std::get<2>(__y));
	}
      else
	{
	  if (__y.index() == 1)
	    return solution_t<_Real>(std::get<2>(__x) * std::get<1>(__y));
	  else
	    return solution_t<_Real>(std::get<2>(__x) * std::get<2>(__y));
	}
    }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator*(const solution_t<_Real>& __x, _Real __y)
    { return operator*(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator*(_Real __x, const solution_t<_Real>& __y)
    { return operator*(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator*(const solution_t<_Real>& __x, std::complex<_Real>& __y)
    { return operator*(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator*(std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator*(solution_t<_Real>(__x), __y); }

  /**
   * division operators...
   */
  template<typename _Real>
    constexpr solution_t<_Real>
    operator/(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      if (__x.index() == 0)
	return __x;
      if (__y.index() == 0)
	return __y;
      else if (__x.index() == 1)
	{
	  if (__y.index() == 1)
	    return solution_t<_Real>(std::get<1>(__x) / std::get<1>(__y));
	  else
	    return solution_t<_Real>(std::get<1>(__x) / std::get<2>(__y));
	}
      else
	{
	  if (__y.index() == 1)
	    return solution_t<_Real>(std::get<2>(__x) / std::get<1>(__y));
	  else
	    return solution_t<_Real>(std::get<2>(__x) / std::get<2>(__y));
	}
    }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator/(const solution_t<_Real>& __x, _Real __y)
    { return operator/(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator/(_Real __x, const solution_t<_Real>& __y)
    { return operator/(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator/(const solution_t<_Real>& __x, std::complex<_Real>& __y)
    { return operator/(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr solution_t<_Real>
    operator/(std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator/(solution_t<_Real>(__x), __y); }

  /**
   * Test for equality and inequality.
   */
  template<typename _Real>
    constexpr bool
    operator==(const solution_t<_Real>& __x,
	       const solution_t<_Real>& __y)
    {
      if (__x.index() == 0 || __y.index() == 0)
	return false;
      else if (__x.index() == __y.index())
	{
	  if (__x.index() == 1)
	    return std::get<1>(__x) == std::get<1>(__y);
	  else
	    return std::get<2>(__x) == std::get<2>(__y);
	}
      else
	return false;
    }

  template<typename _Real>
    bool
    operator==(const solution_t<_Real>& __x, _Real __y)
    { return __x == solution_t<_Real>(__y); }

  template<typename _Real>
    constexpr bool
    operator==(_Real __x, const solution_t<_Real>& __y)
    { return solution_t<_Real>(__x) == __y; }

  template<typename _Real>
    constexpr bool
    operator==(const solution_t<_Real>& __x, const std::complex<_Real>& __y)
    { return __x == solution_t<_Real>(__y); }

  template<typename _Real>
    constexpr bool
    operator==(const std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return solution_t<_Real>(__x) == __y; }

  template<typename _Real>
    constexpr bool
    operator!=(const solution_t<_Real>& __x,
	       const solution_t<_Real>& __y)
    { return !(__x == __y); }

  template<typename _Real>
    bool
    operator!=(const solution_t<_Real>& __x, _Real __y)
    { return !(__x == __y); }

  template<typename _Real>
    constexpr bool
    operator!=(_Real __x, const solution_t<_Real>& __y)
    { return !(__x == __y); }

  template<typename _Real>
    constexpr bool
    operator!=(const solution_t<_Real>& __x, const std::complex<_Real>& __y)
    { return !(__x == __y); }

  template<typename _Real>
    constexpr bool
    operator!=(const std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return !(__x == __y); }

  /**
   * Lexicographic order of solutions as complex numbers.
   * Null solutions compare as less than except to another null solution.
   *
   * A tribool might be a good thing for this when either
   * of the solutions is null.
   */
  template<typename _Real>
    constexpr bool
    operator<(const solution_t<_Real>& __x, const solution_t<_Real>& __y)
    {
      if (__x.index() == 0 && __y.index() == 0)
	return false;
      else if (__x.index() == 0)
	return true;
      else if (__y.index() == 0)
	return false;
      else
	{
	  const auto __rex = real(__x);
	  const auto __rey = real(__y);
	  if (__rex < __rey)
	    return true;
	  else if (__rex == __rey)
	    return imag(__x) < imag(__y);
	  else
	    return false;
	}
    }

  template<typename _Real>
    constexpr bool
    operator<(const solution_t<_Real>& __x, _Real __y)
    { return operator<(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr bool
    operator<(_Real __x, const solution_t<_Real>& __y)
    { return operator<(solution_t<_Real>(__x), __y); }

  template<typename _Real>
    constexpr bool
    operator<(const solution_t<_Real>& __x, const std::complex<_Real>& __y)
    { return operator<(__x, solution_t<_Real>(__y)); }

  template<typename _Real>
    constexpr bool
    operator<(const std::complex<_Real>& __x, const solution_t<_Real>& __y)
    { return operator<(solution_t<_Real>(__x), __y); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

  /**
   * Output a solution to a stream.
   */
  template<typename _Char, typename _Real>
    std::basic_ostream<_Char>&
    operator<<(std::basic_ostream<_Char>& __out,
	       const __gnu_cxx::solution_t<_Real>& __sln)
    {
      const auto __idx = __sln.index();
      if (__idx == 0)
	__out << "null";
      else if (__idx == 1)
	__out << std::get<1>(__sln);
      else if (__idx == 2)
	__out << std::get<2>(__sln);
      return __out;
    }

#endif // C++17

#endif // SOLUTION_H
