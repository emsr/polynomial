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
 * @file solver_low_degree.h
 *
 * This file is a GNU extension to the Standard C++ Library.
 *
 * This file contains the definitions of free functions for solving
 * quadratic, cubic, and quartic equations with real coefficients.
 */


/**
 * @def  _EXT_SOLVER_LOW_DEGREE_TCC
 *
 * @brief  A guard for the low-degree polynomial solver functions header.
 */
#ifndef _EXT_SOLVER_LOW_DEGREE_TCC
#define _EXT_SOLVER_LOW_DEGREE_TCC 1

#pragma GCC system_header

#if __cplusplus < 201402L
# include <bits/c++0x_warning.h>
#else

#include <ext/math_const.h>

namespace __gnu_cxx //_GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   * Refine a solution using the Newton method:
   * @f[
   *    x_{n+1} - x_n = -\frac{P(x_n)}{P'(x_n)}
   * @f]
   */
  template<std::size_t _Dim, typename _Iter, typename _NumTp>
    _NumTp
    __refine_solution_newton(_NumTp __z, const _Iter& _CC)
    {
      for (int __i = 0; __i < 3; ++__i)
	{
	  auto __f = _NumTp(_CC[_Dim - 1]);
	  for (std::size_t __j = _Dim - 1; __j > 0; --__j)
	    __f = _NumTp(_CC[__j - 1]) + __z * __f;

	  auto __df = _NumTp((_Dim - 1) * _CC[_Dim - 1]);
	  for (std::size_t __j = _Dim - 1; __j > 1; --__j)
	    __df = _NumTp((__j - 1) * _CC[__j - 1]) + __z * __df;

	  const auto __del = __f / __df;
	  __z -= __del;
	}
      return __z;
    }

  /**
   * Refine a solution using the Halley method:
   * @f[
   *    x_{n+1} - x_n = -\frac{2 P(x_n) P'(x_n)}
   *                    {2 [P'(x_n)]^2 - P(x_n) P''(x_n)}
   * @f]
   * This form indicates the close relationship to the Newton method:
   * @f[
   *    x_{n+1} - x_n = -\frac{P'(x_n)}
   *                    {P'(x_n) - [P(x_n) P''(x_n)]/[2P'(x_n)]}
   * @f]
   */
  template<std::size_t _Dim, typename _Iter, typename _NumTp>
    _NumTp
    __refine_solution_halley(_NumTp __z, const _Iter& _CC)
    {
      for (int __i = 0; __i < 3; ++__i)
	{
	  auto __f = _NumTp(_CC[_Dim - 1]);
	  for (std::size_t __j = _Dim - 1; __j > 0; --__j)
	    __f = _NumTp(_CC[__j - 1]) + __z * __f;

	  auto __df = _NumTp((_Dim - 1) * _CC[_Dim - 1]);
	  for (std::size_t __j = _Dim - 1; __j > 1; --__j)
	    __df = _NumTp((__j - 1) * _CC[__j - 1]) + __z * __df;

	  auto __d2f = _NumTp((_Dim - 2) * (_Dim - 1) * _CC[_Dim - 1]);
	  for (std::size_t __j = _Dim - 1; __j > 2; --__j)
	    __d2f = _NumTp((__j - 2) * (__j - 1) * _CC[__j - 1]) + __z * __d2f;

	  const auto __del = _NumTp{2} * __f * __df
			   / (_NumTp{2} * __df * __df - __f * __d2f);

	  __z -= __del;
	}
      return __z;
    }

  template<std::size_t _Dim, typename _Iter, typename _Real>
    void
    __refine_solutions(std::array<solution_t<_Real>, _Dim - 1>& _ZZ, const _Iter& _CC)
    {
      for (std::size_t __i = 0; __i < _Dim - 1; ++__i)
	{
	  if (_ZZ[__i].index() == 0)
	    continue;
	  else if (_ZZ[__i].index() == 1)
	    _ZZ[__i] = __refine_solution_newton<_Dim>(std::get<1>(_ZZ[__i]), _CC);
	  else if (_ZZ[__i].index() == 2)
	    _ZZ[__i] = __refine_solution_newton<_Dim>(std::get<2>(_ZZ[__i]), _CC);
	}
    }

  /**
   * @brief Finds the roots of a quadratic equation of the form:
   * @f[
   *    a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   * for real coefficients @f$ a_k @f$.
   *
   * For non-degenerate coefficients two roots are returned:
   * Either the roots are real or the roots are a complex conjugate pair.
   *
   * If the quadratic coefficient @f$ a_2 @f$ is zero (degenerate case)
   * at most one valid root is returned.
   * If the linear coefficient @f$ a_1 @f$ is also zero
   * no valid root is returned.
   *
   * @param[in] _CC Array that contains the three coefficients
   *                  of the quadratic equation.
   */
  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 2>
    __quadratic(const _Iter& _CC)
    {
      std::array<solution_t<_Real>, 2> _ZZ;

      if (_CC[2] == _Real{0})
	{
	  // Equation is linear (or completely degenerate).
	  if (_CC[1] == _Real{0})
	    return _ZZ;
	  else
	    {
	      _ZZ[0] = -_CC[0] / _CC[1];
	      return _ZZ;
	    }
	}
      else if (_CC[0] == _Real{0})
	{
	  _ZZ[0] = _Real{0};
	  if (_CC[2] == _Real{0})
	    return _ZZ;
	  else
	    {
	      _ZZ[1] = -_CC[1] / _CC[2];
	      return _ZZ;
	    }
	}
      else
	{
	  // The discriminant of a quadratic equation
	  const auto _QQ = _CC[1] * _CC[1] - _Real{4} * _CC[2] * _CC[0];

	  if (_QQ < _Real{0})
	    {
	      // The roots are complex conjugates.
	      const auto _ReZZ = -_CC[1] / (_Real{2} * _CC[2]);
	      const auto _ImZZ = std::sqrt(std::abs(_QQ)) / (_Real{2} * _CC[2]);
	      _ZZ[0] = std::complex<_Real>(_ReZZ, -_ImZZ);
	      _ZZ[1] = std::complex<_Real>(_ReZZ, _ImZZ);
	    }
	  else
	    {
	      // The roots are real.
	      _Real __temp = -(_CC[1]
			+ std::copysign(std::sqrt(_QQ), _CC[1])) / _Real{2};
	      _ZZ[0] = __temp / _CC[2];
	      _ZZ[1] = _CC[0] / __temp;
	    }
	}

      return _ZZ;
    }


  /**
   * @brief Finds the roots of a cubic equation of the form:
   * @f[
   *    a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   * for real coefficients @f$ a_k @f$.
   *
   * In the non-degenerate case there are three roots:
   * - All three roots are real
   * - One root is real and the other two are a complex conjugate pair
   *
   * If the cubic coefficient @f$ a_3 @f$ is zero (degenerate case)
   * the problem is referred to the quadratic solver to return, at most,
   * two valid roots.
   *
   * @param[in] _CC Array that contains the four coefficients
   *                  of the cubic equation
   */
  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 3>
    __cubic(const _Iter& _CC)
    {
      using std::experimental::make_array;

      std::array<solution_t<_Real>, 3> _ZZ;

      if (_CC[3] == _Real{0})
	{
	  // Last root is null, remaining equation is quadratic.
	  const auto _ZZ2 = __quadratic<_Real>(_CC);
	  _ZZ[0] = _ZZ2[0];
	  _ZZ[1] = _ZZ2[1];
	}
      else if (_CC[0] == _Real{0})
	{
	  // First root is zero, remaining equation is quadratic.
	  _ZZ[0] = _Real{0};
	  const auto _ZZ2 = __quadratic<_Real>(make_array(_CC[1], _CC[2],
							  _CC[3]));
	  _ZZ[1] = _ZZ2[0];
	  _ZZ[2] = _ZZ2[1];
	}
      else
	{
	  // Normalize cubic equation coefficients.
	  std::array<_Real, 4> _AA3;
	  _AA3[3] = _Real{1};
	  _AA3[2] = _CC[2] / _CC[3];
	  _AA3[1] = _CC[1] / _CC[3];
	  _AA3[0] = _CC[0] / _CC[3];

	  const auto _S_2pi = _Real{2} * __gnu_cxx::__const_pi(_CC[0]);
	  const auto _PP = _AA3[2] / _Real{3};
	  const auto _QQ = (_AA3[2] * _AA3[2] - _Real{3} * _AA3[1])
			 / _Real{9};
	  const auto _QQp3 = _QQ * _QQ * _QQ;
	  const auto _RR = (_Real{2} * _AA3[2] * _AA3[2] * _AA3[2]
			  - _Real{9} * _AA3[2] * _AA3[1]
			  + _Real{27} * _AA3[0]) / _Real{54};
	  const auto _RRp2 = _RR * _RR;

	  if (_QQp3 - _RRp2 > _Real{0})
	    {
	      // Calculate the three real roots.
	      const auto __phi = std::acos(_RR / std::sqrt(_QQp3));
	      const auto __fact = -_Real{2} * std::sqrt(_QQ);
	      for (int __i = 0; __i < 3; ++__i)
		_ZZ[__i] = __fact * std::cos((__phi + __i * _S_2pi) / _Real{3}) - _PP;
	    }
	  else
	    {
	      // Calculate the single real root.
	      const auto __fact = std::cbrt(std::abs(_RR)
					  + std::sqrt(_RRp2 - _QQp3));
	      const auto _BB = -std::copysign(__fact + _QQ / __fact, _RR);
	      _ZZ[0] = _BB - _PP;

	      // Find the other two roots which are complex conjugates.
	      std::array<_Real, 3> _AA2;
	      _AA2[2] = _Real{1};
	      _AA2[1] = _BB;
	      _AA2[0] = _BB * _BB - _Real{3} * _QQ;
	      const auto _ZZ2 = __quadratic<_Real>(_AA2);
	      _ZZ[1] = std::get<2>(_ZZ2[0]) - _PP;
	      _ZZ[2] = std::get<2>(_ZZ2[1]) - _PP;
	    }
	}

      return _ZZ;
    }


  /**
   * @brief Finds the roots a quartic equation of the form:
   * @f[
   * 	 a_4 x^4 + a_3 x^3 + a_2 x^2 + a_1 x + a_0 = 0
   * @f]
   * for real coefficients @f$ a_k @f$.
   *
   * In the non-degenerate case there are four roots:
   * - All four roots are real
   * - Two roots real and two complex roots are a complex conjugate pair
   * - Four complex roots in two complex conjugate pairs
   *
   * If the qartic coefficient @f$ a_4 @f$ is zero (degenerate case)
   * the problem is referred to the cubic solver to return, at most,
   * three valid roots.
   *
   * @param[in] _CC Array that contains the five(5) coefficients
   *                  of the quartic equation.
   */
  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 4>
    __quartic(const _Iter& _CC)
    {
      using std::experimental::make_array;

      std::array<solution_t<_Real>, 4> _ZZ;

      if (_CC[4] == _Real{0})
	{
	  const auto _ZZ3 = __cubic<_Real>(_CC);
	  _ZZ[0] = _ZZ3[0];
	  _ZZ[1] = _ZZ3[1];
	  _ZZ[2] = _ZZ3[2];
	}
      else if (_CC[0] == _Real{0})
	{
	  _ZZ[0] = _Real{0};
	  const auto _ZZ3 = __cubic<_Real>(make_array(_CC[1], _CC[2],
						      _CC[3], _CC[4]));
	  _ZZ[1] = _ZZ3[0];
	  _ZZ[2] = _ZZ3[1];
	}
      else if (_CC[3] == _Real{0} && _CC[1] == _Real{0})
	{
	  // Solve the biquadratic equation.
	  std::array<_Real, 3> _AA2{{_CC[0], _CC[2], _CC[4]}};
	  const auto _ZZ2 = __quadratic<_Real>(_AA2);
	  auto __sqrt = [](solution_t<_Real> __z) -> solution_t<_Real>
			{
			  const auto __idx = __z.index();
			  if (__idx == 0)
			    return __z;
			  else if (__idx == 1)
			    {
			      auto __zz = std::get<1>(__z);
			      return __zz < _Real{0}
				   ? solution_t<_Real>(std::sqrt(std::complex<_Real>(__zz)))
				   : solution_t<_Real>(std::sqrt(__zz));
			    }
			  else
			    return solution_t<_Real>(std::sqrt(std::get<2>(__z)));
			};
	  _ZZ[0] = __sqrt(_ZZ2[0]);
	  _ZZ[1] = __sqrt(_ZZ2[1]);
	  _ZZ[2] = -_ZZ[0];
	  _ZZ[3] = -_ZZ[1];
	}
      else
	{
	  // Normalize quartic equation coefficients.
	  std::array<_Real, 5> _AA4;
	  _AA4[4] = _Real{1};
	  _AA4[3] = _CC[3] / _CC[4];
	  _AA4[2] = _CC[2] / _CC[4];
	  _AA4[1] = _CC[1] / _CC[4];
	  _AA4[0] = _CC[0] / _CC[4];

	  // Calculate the coefficients of the resolvent cubic equation.
	  std::array<_Real, 4> _AA3;
	  _AA3[3] = _Real{1};
	  _AA3[2] = -_AA4[2];
	  _AA3[1] = _AA4[3] * _AA4[1] - _Real{4} * _AA4[0];
	  _AA3[0] = _AA4[0] * (_Real{4} * _AA4[2] - _AA4[3] * _AA4[3])
		  - _AA4[1] * _AA4[1];

	  // Find the algebraically largest real root of the cubic equation
	  // Note: A cubic equation has either three real roots or one
	  //       real root and two complex roots that are complex
	  //       conjugates. If there is only a single real root then
	  //       subroutine cubic always returns that single real root
	  //       (and therefore the algebraically largest real root of
	  //       the cubic equation) as root[0].
	  _Real _Z3max;
	  auto _ZZ3 = __cubic<_Real>(_AA3);
	  if (_ZZ3[1].index() == 1 && _ZZ3[2].index() == 1)
            {
	      // There is some horrible bug with swap and this variant.
	      // They may need to hold the same type.
	      if (_ZZ3[0] < _ZZ3[1])
		//std::swap(_ZZ3[0], _ZZ3[1]);
		{
		  const auto __tmp = _ZZ3[0];
		  _ZZ3[0] = _ZZ3[1];
		  _ZZ3[1] = __tmp;
		}
	      if (_ZZ3[0] < _ZZ3[2])
		//std::swap(_ZZ3[0], _ZZ3[2]);
		{
		  const auto __tmp = _ZZ3[0];
		  _ZZ3[0] = _ZZ3[2];
		  _ZZ3[2] = __tmp;
		}
	      _Z3max = std::get<1>(_ZZ3[0]);
            }
	  else
	    _Z3max = std::get<1>(_ZZ3[0]);

	  // Calculate the coefficients for the two quadratic equations
	  const auto __capa = _Real{0.5L} * _AA4[3];
	  const auto __capb = _Real{0.5L} * _Z3max;
	  const auto __capc = std::sqrt(__capa * __capa - _AA4[2] + _Z3max);
	  const auto __capd = std::sqrt(__capb * __capb - _AA4[0]);
	  const auto __cp = __capa + __capc;
	  const auto __cm = __capa - __capc;
	  auto __dp = __capb + __capd;
	  auto __dm = __capb - __capd;
	  const auto __t1 = __cp * __dm + __cm * __dp;
	  const auto __t2 = __cp * __dp + __cm * __dm;
	  if (std::abs(__t2 - _AA4[1]) < std::abs(__t1 - _AA4[1]))
	    std::swap(__dp, __dm);

	  // Coefficients for the first quadratic equation and find the roots.
	  std::array<_Real, 3> _AA2;
	  _AA2[2] = _Real{1};
	  _AA2[1] = __cp;
	  _AA2[0] = __dp;
	  const auto _ZZ2p = __quadratic<_Real>(_AA2);
	  _ZZ[0] = _ZZ2p[0];
	  _ZZ[1] = _ZZ2p[1];

	  // Coefficients for the second quadratic equation and find the roots.
	  _AA2[2] = _Real{1};
	  _AA2[1] = __cm;
	  _AA2[0] = __dm;
	  const auto _ZZ2m = __quadratic<_Real>(_AA2);
	  _ZZ[2] = _ZZ2m[0];
	  _ZZ[3] = _ZZ2m[1];
	}

      return _ZZ;
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // C++14

#endif // _EXT_SOLVER_LOW_DEGREE_TCC
