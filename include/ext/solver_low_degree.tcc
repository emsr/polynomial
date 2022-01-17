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
 * @file solver_low_degree.h
 *
 * This file is a GNU extension to the Standard C++ Library.
 *
 * This file contains the definitions of free functions for solving
 * quadratic, cubic, and quartic equations with real coefficients.
 */


/**
 * @def  SOLVER_LOW_DEGREE_TCC
 *
 * @brief  A guard for the low-degree polynomial solver functions header.
 */
#ifndef SOLVER_LOW_DEGREE_TCC
#define SOLVER_LOW_DEGREE_TCC 1

namespace emsr
{

  /**
   * Refine a solution using the Newton method:
   * @f[
   *    x_{n+1} - x_n = -\frac{P(x_n)}{P'(x_n)}
   * @f]
   */
  template<std::size_t _Dim, typename _Iter, typename _NumTp>
    _NumTp
    refine_solution_newton(_NumTp z, const _Iter& _CC)
    {
      for (int i = 0; i < 3; ++i)
	{
	  auto f = _NumTp(_CC[_Dim - 1]);
	  for (std::size_t j = _Dim - 1; j > 0; --j)
	    f = _NumTp(_CC[j - 1]) + z * f;

	  auto df = _NumTp((_Dim - 1) * _CC[_Dim - 1]);
	  for (std::size_t j = _Dim - 1; j > 1; --j)
	    df = _NumTp((j - 1) * _CC[j - 1]) + z * df;

	  const auto del = f / df;
	  z -= del;
	}
      return z;
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
    refine_solution_halley(_NumTp z, const _Iter& _CC)
    {
      for (int i = 0; i < 3; ++i)
	{
	  auto f = _NumTp(_CC[_Dim - 1]);
	  for (std::size_t j = _Dim - 1; j > 0; --j)
	    f = _NumTp(_CC[j - 1]) + z * f;

	  auto df = _NumTp((_Dim - 1) * _CC[_Dim - 1]);
	  for (std::size_t j = _Dim - 1; j > 1; --j)
	    df = _NumTp((j - 1) * _CC[j - 1]) + z * df;

	  auto d2f = _NumTp((_Dim - 2) * (_Dim - 1) * _CC[_Dim - 1]);
	  for (std::size_t j = _Dim - 1; j > 2; --j)
	    d2f = _NumTp((j - 2) * (j - 1) * _CC[j - 1]) + z * d2f;

	  const auto del = _NumTp{2} * f * df
			   / (_NumTp{2} * df * df - f * d2f);

	  z -= del;
	}
      return z;
    }

  template<std::size_t _Dim, typename _Iter, typename Real>
    void
    refine_solutions(std::array<solution_t<Real>, _Dim - 1>& _ZZ, const _Iter& _CC)
    {
      for (std::size_t i = 0; i < _Dim - 1; ++i)
	{
	  if (_ZZ[i].index() == 0)
	    continue;
	  else if (_ZZ[i].index() == 1)
	    _ZZ[i] = refine_solution_newton<_Dim>(std::get<1>(_ZZ[i]), _CC);
	  else if (_ZZ[i].index() == 2)
	    _ZZ[i] = refine_solution_newton<_Dim>(std::get<2>(_ZZ[i]), _CC);
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
  template<typename Real, typename _Iter>
    std::array<solution_t<Real>, 2>
    quadratic(const _Iter& _CC)
    {
      std::array<solution_t<Real>, 2> _ZZ;

      if (_CC[2] == Real{0})
	{
	  // Equation is linear (or completely degenerate).
	  if (_CC[1] == Real{0})
	    return _ZZ;
	  else
	    {
	      _ZZ[0] = -_CC[0] / _CC[1];
	      return _ZZ;
	    }
	}
      else if (_CC[0] == Real{0})
	{
	  _ZZ[0] = Real{0};
	  if (_CC[2] == Real{0})
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
	  const auto _QQ = _CC[1] * _CC[1] - Real{4} * _CC[2] * _CC[0];

	  if (_QQ < Real{0})
	    {
	      // The roots are complex conjugates.
	      const auto _ReZZ = -_CC[1] / (Real{2} * _CC[2]);
	      const auto _ImZZ = std::sqrt(std::abs(_QQ)) / (Real{2} * _CC[2]);
	      _ZZ[0] = std::complex<Real>(_ReZZ, -_ImZZ);
	      _ZZ[1] = std::complex<Real>(_ReZZ, _ImZZ);
	    }
	  else
	    {
	      // The roots are real.
	      Real temp = -(_CC[1]
			+ std::copysign(std::sqrt(_QQ), _CC[1])) / Real{2};
	      _ZZ[0] = temp / _CC[2];
	      _ZZ[1] = _CC[0] / temp;
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
  template<typename Real, typename Iter>
    std::array<solution_t<Real>, 3>
    cubic(const Iter& CC)
    {
      std::array<solution_t<Real>, 3> ZZ;

      if (CC[3] == Real{0})
	{
	  // Last root is null, remaining equation is quadratic.
	  const auto ZZ2 = quadratic<Real>(CC);
	  ZZ[0] = ZZ2[0];
	  ZZ[1] = ZZ2[1];
	}
      else if (CC[0] == Real{0})
	{
	  // First root is zero, remaining equation is quadratic.
	  ZZ[0] = Real{0};
	  const auto ZZ2 = quadratic<Real>(CC[1], CC[2], CC[3]);
	  ZZ[1] = ZZ2[0];
	  ZZ[2] = ZZ2[1];
	}
      else
	{
	  // Normalize cubic equation coefficients.
	  std::array<Real, 4> AA3;
	  AA3[3] = Real{1};
	  AA3[2] = CC[2] / CC[3];
	  AA3[1] = CC[1] / CC[3];
	  AA3[0] = CC[0] / CC[3];

	  const auto S_2pi
	    = 2 * Real{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
	  const auto PP = AA3[2] / Real{3};
	  const auto QQ = (AA3[2] * AA3[2] - Real{3} * AA3[1])
			 / Real{9};
	  const auto QQp3 = QQ * QQ * QQ;
	  const auto RR = (Real{2} * AA3[2] * AA3[2] * AA3[2]
			  - Real{9} * AA3[2] * AA3[1]
			  + Real{27} * AA3[0]) / Real{54};
	  const auto RRp2 = RR * RR;

	  if (QQp3 - RRp2 > Real{0})
	    {
	      // Calculate the three real roots.
	      const auto phi = std::acos(RR / std::sqrt(QQp3));
	      const auto fact = -Real{2} * std::sqrt(QQ);
	      for (int i = 0; i < 3; ++i)
		ZZ[i] = fact * std::cos((phi + i * S_2pi) / Real{3}) - PP;
	    }
	  else
	    {
	      // Calculate the single real root.
	      const auto fact = std::cbrt(std::abs(RR)
					  + std::sqrt(RRp2 - QQp3));
	      const auto BB = -std::copysign(fact + QQ / fact, RR);
	      ZZ[0] = BB - PP;

	      // Find the other two roots which are complex conjugates.
	      std::array<Real, 3> AA2;
	      AA2[2] = Real{1};
	      AA2[1] = BB;
	      AA2[0] = BB * BB - Real{3} * QQ;
	      const auto ZZ2 = quadratic<Real>(AA2);
	      ZZ[1] = std::get<2>(ZZ2[0]) - PP;
	      ZZ[2] = std::get<2>(ZZ2[1]) - PP;
	    }
	}

      return ZZ;
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
   * @param[in] CC Array that contains the five(5) coefficients
   *                  of the quartic equation.
   */
  template<typename Real, typename Iter>
    std::array<solution_t<Real>, 4>
    quartic(const Iter& CC)
    {
      std::array<solution_t<Real>, 4> ZZ;

      if (CC[4] == Real{0})
	{
	  const auto ZZ3 = cubic<Real>(CC);
	  ZZ[0] = ZZ3[0];
	  ZZ[1] = ZZ3[1];
	  ZZ[2] = ZZ3[2];
	}
      else if (CC[0] == Real{0})
	{
	  ZZ[0] = Real{0};
	  const auto ZZ3 = cubic<Real>(CC[1], CC[2], CC[3], CC[4]);
	  ZZ[1] = ZZ3[0];
	  ZZ[2] = ZZ3[1];
	}
      else if (CC[3] == Real{0} && CC[1] == Real{0})
	{
	  // Solve the biquadratic equation.
	  std::array<Real, 3> AA2{{CC[0], CC[2], CC[4]}};
	  const auto ZZ2 = quadratic<Real>(AA2);
	  auto sqrt = [](solution_t<Real> z) -> solution_t<Real>
			{
			  const auto idx = z.index();
			  if (idx == 0)
			    return z;
			  else if (idx == 1)
			    {
			      auto zz = std::get<1>(z);
			      return zz < Real{0}
				   ? solution_t<Real>(std::sqrt(std::complex<Real>(zz)))
				   : solution_t<Real>(std::sqrt(zz));
			    }
			  else
			    return solution_t<Real>(std::sqrt(std::get<2>(z)));
			};
	  ZZ[0] = sqrt(ZZ2[0]);
	  ZZ[1] = sqrt(ZZ2[1]);
	  ZZ[2] = -ZZ[0];
	  ZZ[3] = -ZZ[1];
	}
      else
	{
	  // Normalize quartic equation coefficients.
	  std::array<Real, 5> AA4;
	  AA4[4] = Real{1};
	  AA4[3] = CC[3] / CC[4];
	  AA4[2] = CC[2] / CC[4];
	  AA4[1] = CC[1] / CC[4];
	  AA4[0] = CC[0] / CC[4];

	  // Calculate the coefficients of the resolvent cubic equation.
	  std::array<Real, 4> AA3;
	  AA3[3] = Real{1};
	  AA3[2] = -AA4[2];
	  AA3[1] = AA4[3] * AA4[1] - Real{4} * AA4[0];
	  AA3[0] = AA4[0] * (Real{4} * AA4[2] - AA4[3] * AA4[3])
		  - AA4[1] * AA4[1];

	  // Find the algebraically largest real root of the cubic equation
	  // Note: A cubic equation has either three real roots or one
	  //       real root and two complex roots that are complex
	  //       conjugates. If there is only a single real root then
	  //       subroutine cubic always returns that single real root
	  //       (and therefore the algebraically largest real root of
	  //       the cubic equation) as root[0].
	  Real Z3max;
	  auto ZZ3 = cubic<Real>(AA3);
	  if (ZZ3[1].index() == 1 && ZZ3[2].index() == 1)
            {
	      // There is some horrible bug with swap and this variant.
	      // They may need to hold the same type.
	      if (ZZ3[0] < ZZ3[1])
		//std::swap(ZZ3[0], ZZ3[1]);
		{
		  const auto tmp = ZZ3[0];
		  ZZ3[0] = ZZ3[1];
		  ZZ3[1] = tmp;
		}
	      if (ZZ3[0] < ZZ3[2])
		//std::swap(ZZ3[0], ZZ3[2]);
		{
		  const auto tmp = ZZ3[0];
		  ZZ3[0] = ZZ3[2];
		  ZZ3[2] = tmp;
		}
	      Z3max = std::get<1>(ZZ3[0]);
            }
	  else
	    Z3max = std::get<1>(ZZ3[0]);

	  // Calculate the coefficients for the two quadratic equations
	  const auto capa = Real{0.5L} * AA4[3];
	  const auto capb = Real{0.5L} * Z3max;
	  const auto capc = std::sqrt(capa * capa - AA4[2] + Z3max);
	  const auto capd = std::sqrt(capb * capb - AA4[0]);
	  const auto cp = capa + capc;
	  const auto cm = capa - capc;
	  auto dp = capb + capd;
	  auto dm = capb - capd;
	  const auto t1 = cp * dm + cm * dp;
	  const auto t2 = cp * dp + cm * dm;
	  if (std::abs(t2 - AA4[1]) < std::abs(t1 - AA4[1]))
	    std::swap(dp, dm);

	  // Coefficients for the first quadratic equation and find the roots.
	  std::array<Real, 3> AA2;
	  AA2[2] = Real{1};
	  AA2[1] = cp;
	  AA2[0] = dp;
	  const auto ZZ2p = quadratic<Real>(AA2);
	  ZZ[0] = ZZ2p[0];
	  ZZ[1] = ZZ2p[1];

	  // Coefficients for the second quadratic equation and find the roots.
	  AA2[2] = Real{1};
	  AA2[1] = cm;
	  AA2[0] = dm;
	  const auto ZZ2m = quadratic<Real>(AA2);
	  ZZ[2] = ZZ2m[0];
	  ZZ[3] = ZZ2m[1];
	}

      return ZZ;
    }

} // namespace emsr

#endif // SOLVER_LOW_DEGREE_TCC
