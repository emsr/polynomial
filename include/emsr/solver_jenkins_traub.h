
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
 * @file solver_jenkins_traub.h Class declaration for the Jenkins-Traub solver.
 */

/**
 * @def  SOLVER_JENKINS_TRAUB_H
 *
 * @brief  A guard for the _JenkinsTraubSolver class header.
 */
#ifndef SOLVER_JENKINS_TRAUB_H
#define SOLVER_JENKINS_TRAUB_H 1

#include <emsr/solution.h> // For solution_t

namespace emsr
{

/**
 * A solver for real-coefficient polynomials due to Jenkins and Traub.
 */
template<typename Real>
  class JenkinsTraubSolver
  {
  public:

    JenkinsTraubSolver(const std::vector<Real>& op);
    JenkinsTraubSolver(std::vector<Real>&& op);

    std::vector<solution_t<Real>> solve();

  private:

    enum NormalizationType
    {
      none,
      divide_by_c,
      divide_by_d,
      near_h_root
    };

    void quadratic(Real a, Real b, Real c,
		   solution_t<Real> &z_small, solution_t<Real> &z_large);
    int fxshfr(int l2);
    int iter_quadratic(Real uu, Real vv);
    int iter_real(Real sss, int& iflag);
    NormalizationType init_next_h_poly();
    void next_h_poly(NormalizationType type);
    std::pair<Real, Real> quadratic_coefficients(NormalizationType type);
    void remquo_quadratic(int n, Real u, Real v,
			  std::vector<Real>& poly, std::vector<Real>& quot,
			  Real& a, Real& b);

    static constexpr auto s_eps = std::numeric_limits<Real>::epsilon();
    static constexpr auto s_base = Real{std::numeric_limits<Real>::radix};
    static constexpr auto s_tiny = s_eps * s_eps * s_eps; // 1.0e-50; //std::numeric_limits<Real>::min();
    static constexpr auto s_huge = std::numeric_limits<Real>::max();
    //static constexpr auto s_low = s_tiny / s_eps;
    const Real s_low = s_tiny / s_eps;
    static constexpr auto s_pi = Real{3.1415'92653'58979'32384'62643'38327'95028'84195e+0L};
    static constexpr auto s_rotation = Real{94} * s_pi / Real{180};

    int m_max_iter_quadratic = 20;
    Real m_min_log_deriv = Real{0.005L};
    int m_max_iter_real = 10;
    // Epsilon parameter.
    Real m_are = s_eps;
    // Epsilon parameter.
    Real m_mre = s_eps;

    std::vector<Real> m_P;
    std::vector<Real> m_P_quot;
    std::vector<Real> m_H, m_H_quot, m_H_save;
    Real m_sr, m_si;
    Real m_u, m_v;
    Real m_a;
    Real m_b;
    Real m_c;
    Real m_d;
    Real m_e;
    Real m_f;
    Real m_g;
    Real m_h;
    Real m_a1;
    Real m_a2;
    Real m_a3;
    Real m_a7;
    solution_t<Real> m_z_small;
    solution_t<Real> m_z_large;
    int m_order;
    bool m_zerok;
    int m_num_iters = 0;
  };

} // namespace emsr

#include <emsr/solver_jenkins_traub.tcc>
#include <emsr/solver_jenkins_traub_complex.tcc>

#endif // SOLVER_JENKINS_TRAUB_H
