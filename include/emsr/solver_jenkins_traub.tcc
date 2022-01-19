
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
 * @file solver_jenkins_traub.tcc Class declaration for
 * the Jenkins-Traub solver.
 */

/**
 * @def  SOLVER_JENKINS_TRAUB_TCC
 *
 * @brief  A guard for the _JenkinsTraubSolver class header.
 */
#ifndef SOLVER_JENKINS_TRAUB_TCC
#define SOLVER_JENKINS_TRAUB_TCC 1

#include <emsr/solver_low_degree.h>

namespace emsr
{

/**
 * Constructor from input polynomial.
 */
template<typename Real>
  JenkinsTraubSolver<Real>::
  JenkinsTraubSolver(const std::vector<Real>& op)
  : m_P(op)
  {
    if (this->m_P.size() == 0)
      throw std::domain_error("Polynomial degree must be at least 1.");

    // Algorithm fails if the leading coefficient is zero.
    // We could erase leading-order zero coefficients.
    if (this->m_P[0] == Real{0})
      throw std::domain_error("Leading coefficient must be nonzero.");

    const auto degree = this->m_P.size() - 1;
    this->m_order = degree;
    this->m_P_quot.resize(degree + 1);
    this->m_H.resize(degree + 1);
    this->m_H_quot.resize(degree + 1);
    this->m_H_save.resize(degree + 1);
  }

/**
 *
 */
template<typename Real>
  std::vector<Solution<Real>>
  JenkinsTraubSolver<Real>::solve()
  {
    // Initialization of constants for shift rotation.
    auto xx = 1 / Real{1.4142'13562'37309'50488'01688'72420'96980'78569e+0L};
    auto yy = -xx;
    const auto cosr = std::cos(s_rotation);
    const auto sinr = std::sin(s_rotation);

    std::vector<Solution<Real>> zero;
    zero.reserve(this->m_P.size());

    // Remove the zeros at the origin, if any.
    while (this->m_P[this->m_order] == Real{0})
      {
	zero.push_back(Real{0});
	--this->m_order;
      }
    if (this->m_order < 1)
      return zero;

    std::vector<Real> pt(this->m_order + 1);
    std::vector<Real> _H_temp(this->m_order + 1);

    while (true)
      {
	// Start the algorithm for one zero.
	this->m_num_iters = 0;
	if (this->m_order == 1)
	  {
	    zero.push_back(-this->m_P[1] / this->m_P[0]);
	    --this->m_order;
	    return zero;
	  }
	// Calculate the final zero or pair of zeros.
	if (this->m_order == 2)
	  {
	    Solution<Real> z_small, z_large;
	    this->quadratic(this->m_P[0], this->m_P[1], this->m_P[2],
			    z_small, z_large);
	    if (z_small.index() != 0)
	      {
		zero.push_back(z_small);
		--this->m_order;
	      }
	    if (z_large.index() != 0)
	      {
		zero.push_back(z_large);
		--this->m_order;
	      }
	    return zero;
	  }

	// Find largest and smallest moduli of coefficients.
	auto a_max = Real{0};
	auto a_min = s_huge;
	for (int i = 0; i <= this->m_order; ++i)
	  {
	    const auto x = std::abs(this->m_P[i]);
	    if (x > a_max)
	      a_max = x;
	    if (x != Real{0} && x < a_min)
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
	if (scale > Real{1} && s_huge / scale < a_max)
	  rescale = false;
	if (scale <= Real{1} && a_max < Real{10})
	  rescale = false;

	if (rescale)
	  {
	    // Scale polynomial.
	    if (scale == Real{0})
	      scale = s_tiny;
	    const auto l = std::ilogb(scale);
	    const auto factor = std::pow(s_base, l);
	    if (factor != Real{1})
	      for (int i = 0; i <= this->m_order; ++i)
		this->m_P[i] *= factor;
	  }

	// Compute lower bound on moduli of roots.
	for (int i = 0; i <= this->m_order; ++i)
	  pt[i] = std::abs(this->m_P[i]);
	pt[this->m_order] = -pt[this->m_order];
	// Compute upper estimate of bound.
	auto x = std::exp((std::log(-pt[this->m_order])
			    - std::log(pt[0])) / Real(this->m_order));
	// If Newton step at the origin is better, use it.	
	if (pt[this->m_order - 1] != Real{0})
	  {
	    const auto xm = -pt[this->m_order] / pt[this->m_order - 1];
	    if (xm < x)
	      x = xm;
	  }
	// Chop the interval (0,x) until ff <= 0.
	while (true)
	  {
	    auto xm = x * Real{0.1L};
	    auto ff = pt[0];
	    for (int i = 1; i <= this->m_order; ++i)
	      ff = ff * xm + pt[i];
	    if (ff <= Real{0})
	      break;
	    x = xm;
	  }
	// Do Newton interation until x converges to two decimal places.
	auto dx = x;
	while (std::abs(dx / x) > this->m_min_log_deriv)
	  {
	    auto ff = pt[0];
	    auto df = ff;
	    for (int i = 1; i < this->m_order; ++i)
	      {
		ff = ff * x + pt[i];
		df = df * x + ff;
	      }
	    ff = ff * x + pt[this->m_order];
	    dx = ff / df;
	    x -= dx;
	    ++this->m_num_iters;
	  }
	const auto bound = x;
	// Compute the derivative as the initial _H polynomial
	// and do 5 steps with no shift.
	const auto nm1 = this->m_order - 1;
	for (int i = 1; i < this->m_order; ++i)
	  this->m_H[i] = Real(this->m_order - i) * this->m_P[i]
			/ Real(this->m_order);
	this->m_H[0] = this->m_P[0];
	const auto aa = this->m_P[this->m_order];
	const auto bb = this->m_P[this->m_order - 1];
	this->m_zerok = (this->m_H[this->m_order - 1] == Real{0});
	for(int jj = 0; jj < 5; ++jj)
	  {
	    ++this->m_num_iters;
	    auto cc = this->m_H[this->m_order - 1];
	    if (!this->m_zerok)
	      {
		// Use a scaled form of recurrence if value of H at 0
		// is nonzero.	
		const auto t = -aa / cc;
		for (int i = 0; i < nm1; ++i)
		  {
		    const auto j = this->m_order - i - 1;
		    this->m_H[j] = t * this->m_H[j - 1] + this->m_P[j];
		  }
		this->m_H[0] = this->m_P[0];
		this->m_zerok = (std::abs(this->m_H[this->m_order - 1])
			    <= Real{10} * s_eps * std::abs(bb));
	    }
	    else
	      {
		// Use unscaled form of recurrence.
		for (int i = 0; i < nm1; ++i)
		  {
		    const auto j = this->m_order - i - 1;
		    this->m_H[j] = this->m_H[j - 1];
		  }
		this->m_H[0] = Real{0};
		this->m_zerok = (this->m_H[this->m_order - 1] == Real{0});
	      }
	  }
	// Save H for restarts with new shifts.
	_H_temp = this->m_H;

	// Loop to select the quadratic corresponding to each new shift.
	for (int count = 0; count < 20; ++count)
	  {
	    /*  Quadratic corresponds to a Real shift to a	
	     *  non-real point and its complex conjugate. The point
	     *  has modulus bound and amplitude rotated by 94 degrees
	     *  from the previous shift.
	     */
	    const auto xxx = cosr * xx - sinr * yy;
	    yy = sinr * xx + cosr * yy;
	    auto xx = xxx;
	    this->m_sr = bound * xx;
	    this->m_si = bound * yy;
	    this->m_u = -Real{2} * this->m_sr;
	    this->m_v = bound;
	    auto num_zeros = this->fxshfr(20 * (count + 1));
	    bool cycle = false;
	    if (num_zeros != 0)
	      {
	      /*  The second stage jumps directly to one of the third
	       *  stage iterations and returns here if successful.
	       *  Deflate the polynomial, store the zero or zeros
	       *  and return to the main algorithm.
	       */
		zero.push_back(this->m_z_small);
		this->m_order -= num_zeros;
		this->m_P = this->m_P_quot;
		if (num_zeros != 1)
		  zero.push_back(this->m_z_large);
		cycle = true;
		break;
	      }
	    if (cycle)
	      continue;

	    // If the iteration is unsuccessful another quadratic
	    // is chosen after restoring H.
	    this->m_H = _H_temp;
	 }
      }
  }


/**
 * Computes up to L2 fixed shift H-polynomials, testing for convergence
 * in the linear or quadratic case.
 * Initiates one of the variable shift iterations and returns
 * the number of zeros found.
 */
template<typename Real>
  int
  JenkinsTraubSolver<Real>::fxshfr(int l2)
  {
    Real ts_old{}, tv_old{};
    int iflag;

    int num_zeros = 0;

    auto betav = Real{0.25};
    auto betas = Real{0.25};
    auto ss_old = this->m_sr;
    auto vv_old = this->m_v;
    // Evaluate polynomial by synthetic division.
    this->remquo_quadratic(this->m_order, this->m_u, this->m_v,
			   this->m_P, this->m_P_quot,
			   this->m_a, this->m_b);
    auto type = this->init_next_h_poly();
    for (int j = 0; j < l2; ++j)
      {
	// Calculate next H polynomial and estimate v.
	this->next_h_poly(type);
	type = this->init_next_h_poly();
	auto [ui, vi] = this->quadratic_coefficients(type);
	auto vv = vi;
	// Estimate s.
	auto ss = Real{0};
	if (this->m_H[this->m_order - 1] != Real{0})
	  ss = -this->m_P[this->m_order] / this->m_H[this->m_order - 1];
	auto tv = Real{1};
	auto ts = Real{1};
	if (j == 0 || type == near_h_root)
	  {
	    vv_old = vv;
	    ss_old = ss;
	    tv_old = tv;
	    ts_old = ts;
	    continue;
	  }

	// Compute relative measures of convergence of s and v sequences.
	if (vv != Real{0})
	  tv = std::abs((vv - vv_old) / vv);
	if (ss != Real{0})
	  ts = std::abs((ss - ss_old) / ss);

	// If decreasing, multiply two most recent convergence measures.
	const auto tvv = tv < tv_old ? tv * tv_old : Real{1};
	const auto tss = ts < ts_old ? ts * ts_old : Real{1};

	// Compare with convergence criteria.
	const auto vpass = tvv < betav;
	const auto spass = tss < betas;
	if (!(spass || vpass))
	  {
	    vv_old = vv;
	    ss_old = ss;
	    tv_old = tv;
	    ts_old = ts;
	    continue;
	  }

	// At least one sequence has passed the convergence test.
	// Store variables before iterating.
	const auto u_save = this->m_u;
	const auto v_save = this->m_v;
	this->m_H_save = this->m_H;
	const auto s = ss;

	// Choose iteration according to the fastest converging sequence.
	auto vtry = false;
	auto stry = false;
	if ((spass && !vpass) || tss < tvv)
	  goto TRY_LINEAR;

  TRY_QUADRATIC:
	num_zeros = this->iter_quadratic(ui, vi);
	if (num_zeros > 0)
	  return num_zeros;
	// Quadratic iteration has failed. Flag that it has
	// been tried and decrease the convergence criterion.
	vtry = true;
	betav *= Real{0.25};
	// Try linear iteration if it has not been tried and
	// the S sequence is converging.
	if (stry || !spass)
	  goto RESTORE_VARS;
	this->m_H = this->m_H_save;

  TRY_LINEAR:
	num_zeros = this->iter_real(s, iflag);
	if (num_zeros > 0)
	  return num_zeros;
	// Linear iteration has failed. Flag that it has been
	// tried and decrease the convergence criterion.
	stry = true;
	betas *= Real{0.25};
	if (iflag == 0)
	  goto RESTORE_VARS;
	// If linear iteration signals an almost real
	// zero attempt quadratic iteration.
	ui = -Real{2} * s;
	vi = s * s;
	goto TRY_QUADRATIC;

  RESTORE_VARS:
	// Restore variables.
	this->m_u = u_save;
	this->m_v = v_save;
	this->m_H = this->m_H_save;

	// Try quadratic iteration if it has not been tried
	// and the V sequence is converging.
	if (!vtry && vpass)
	  goto TRY_QUADRATIC;

	// Recompute polynomial quotient and remainder
        // to continue the second stage.
	this->remquo_quadratic(this->m_order, this->m_u, this->m_v,
			       this->m_P, this->m_P_quot,
			       this->m_a, this->m_b);
	type = this->init_next_h_poly();
      }
    return num_zeros;
  }


/**
 * Variable-shift H-polynomial iteration for a quadratic factor
 * converges only if the zeros are equimodular or nearly so.
 * @param uu The linear coefficient of the starting quadratic equation.
 * @param vv The constant  coefficient of the starting quadratic equation.
 * @return The number of zeros found.
 */
template<typename Real>
  int
  JenkinsTraubSolver<Real>::iter_quadratic(Real uu, Real vv)
  {
    Real mp, mp_old{}, ee, relstp = std::sqrt(s_eps), t, zm;
    NormalizationType type;

    int num_zeros = 0;
    bool tried = false;
    this->m_u = uu;
    this->m_v = vv;
    int j = 0;

    while (true)
      {
	++this->m_num_iters;
	this->quadratic(Real{1}, this->m_u, this->m_v,
			this->m_z_small, this->m_z_large);
	// Return if roots of the quadratic are real and not
	// close to multiple or nearly equal and of opposite sign.
	if (std::abs(std::abs(real(this->m_z_small))
		   - std::abs(real(this->m_z_large)))
	       > Real{0.01L} * std::abs(real(this->m_z_large)))
	  return num_zeros;
	// Evaluate polynomial by quadratic synthetic division.
	this->remquo_quadratic(this->m_order, this->m_u, this->m_v,
			       this->m_P, this->m_P_quot, this->m_a, this->m_b);
	mp = std::abs(this->m_a - real(this->m_z_small) * this->m_b)
	   + std::abs(imag(this->m_z_small) * this->m_b);
	// Compute a rigorous bound on the rounding error in evaluating _P.
	zm = std::sqrt(std::abs(this->m_v));
	ee = Real{2} * std::abs(this->m_P_quot[0]);
	t = -real(this->m_z_small) * this->m_b;
	for (int i = 1; i < this->m_order; ++i)
	  ee = ee * zm + std::abs(this->m_P_quot[i]);
	ee = ee * zm + std::abs(this->m_a + t);
	ee *= (Real{5} * this->m_mre + Real{4} * this->m_are);
	ee -= (Real{5} * this->m_mre + Real{2} * this->m_are)
	    * (std::abs(this->m_a + t) + std::abs(this->m_b) * zm);
	ee += Real{2} * this->m_are * std::abs(t);
	// Iteration has converged sufficiently if the
	// polynomial value is less than 20 times this bound.
	if (mp <= Real{20} * ee)
	  {
	    num_zeros = 2;
	    return num_zeros;
	  }
	++j;
	// Stop iteration after 20 steps.
	if (j > this->m_max_iter_quadratic)
	  return num_zeros;
	if (j < 2 || relstp > Real{0.01L} || mp < mp_old || tried)
	  {
	    mp_old = mp;
	    // Calculate next H polynomial and new u and v.
	    type = this->init_next_h_poly();
	    this->next_h_poly(type);
	    type = this->init_next_h_poly();
	    const auto [ui, vi] = this->quadratic_coefficients(type);
	    // If vi is zero the iteration is not converging.
	    if (vi == Real{0})
	      return num_zeros;
	    relstp = std::abs((vi - this->m_v) / vi);
	    this->m_u = ui;
	    this->m_v = vi;
	    continue;
	  }
	// A cluster appears to be stalling the convergence.
	// Five fixed shift steps are taken with a u, v close to the cluster.
	if (relstp < s_eps)
	  relstp = s_eps;
	relstp = std::sqrt(relstp);
	this->m_u -= this->m_u * relstp;
	this->m_v += this->m_v * relstp;
	this->remquo_quadratic(this->m_order, this->m_u, this->m_v,
			       this->m_P, this->m_P_quot,
			       this->m_a, this->m_b);
	for (int i = 0; i < 5; ++i)
	  {
	    type = this->init_next_h_poly();
	    this->next_h_poly(type);
	  }
	tried = true;
	j = 0;
      }
  }


/**
 * Variable-shift H polynomial iteration for a real zero.
 * @param sss Starting iterate.
 * @param iflag Flag to indicate a pair of zeros near real axis.
 * @return The number of zeros found.
 */
template<typename Real>
  int
  JenkinsTraubSolver<Real>::iter_real(Real sss, int& iflag)
  {
    auto t = Real{0};
    decltype(std::abs(this->m_P[0])) mp_old{};

    int num_zeros = 0;
    auto s = sss;
    iflag = 0;
    int i_real = 0;

    while (true)
      {
	++this->m_num_iters;
	auto pval = this->m_P[0];
	// Evaluate P at s.
	this->m_P_quot[0] = pval;
	for (int i = 1; i <= this->m_order; ++i)
	  {
	    pval = pval * s + this->m_P[i];
	    this->m_P_quot[i] = pval;
	  }
	auto mp = std::abs(pval);
	// Compute a rigorous bound on the error in evaluating P.
	const auto ms = std::abs(s);
	auto ee = (this->m_mre / (this->m_are + this->m_mre))
		  * std::abs(this->m_P_quot[0]);
	for (int i = 1; i <= this->m_order; ++i)
	  ee = ee * ms + std::abs(this->m_P_quot[i]);
	// Iteration has converged sufficiently if the polynomial
	// value is less than 20 times this bound.
	if (mp <= Real{20}
		 * ((this->m_are + this->m_mre) * ee - this->m_mre * mp))
	  {
	    num_zeros = 1;
	    this->m_z_small = s;
	    return num_zeros;
	  }
	++i_real;
	// Stop iteration after max_iter_real steps.
	if (i_real > this->m_max_iter_real)
	  return num_zeros;
	else if (i_real < 2
	  || std::abs(t) > Real{0.001L} * std::abs(s - t)
	  || mp < mp_old)
	  {
	    // Return if the polynomial value has increased significantly.
	    mp_old = mp;

	    // Compute t, the next polynomial, and the new iterate.
	    auto hval = this->m_H[0];
	    this->m_H_quot[0] = hval;
	    for (int i = 1; i < this->m_order; ++i)
	      {
		hval = hval * s + this->m_H[i];
		this->m_H_quot[i] = hval;
	      }

	    if (std::abs(hval)
		 <= std::abs(this->m_H[this->m_order - 1]) * Real{10} * s_eps)
	      {
		// Use unscaled form.
		this->m_H[0] = Real{0};
		for (int i = 1; i < this->m_order; ++i)
		  this->m_H[i] = this->m_H_quot[i-1];
	      }
	    else
	      {
		// Use the scaled form of the recurrence if the value
		// of H at s is nonzero.
		t = -pval / hval;
		this->m_H[0] = this->m_P_quot[0];
		for (int i = 1; i < this->m_order; ++i)
		  this->m_H[i] = t * this->m_H_quot[i - 1]
				+ this->m_P_quot[i];
	      }

	    hval = this->m_H[0];
	    for (int i = 1; i < this->m_order; ++i)
	      hval = hval * s + this->m_H[i];
	    auto t = Real{0};
	    if (std::abs(hval)
		 > std::abs(this->m_H[this->m_order - 1] * Real{10} * s_eps))
	      t = -pval / hval;
	    s += t;
	  }
	else
	  {
	    // A cluster of zeros near the real axis has been encountered.
	    // Return with iflag set to initiate a quadratic iteration.
	    iflag = 1;
	    sss = s;
	    return num_zeros;
	  }
      }
  }

/**
 * This routine calculates scalar quantities used to compute
 * the next H-polynomial and new estimates of the quadratic coefficients.
 *
 * @return Flag indicating how the calculations are normalized
 *         to avoid overflow.
 */
template<typename Real>
  typename JenkinsTraubSolver<Real>::NormalizationType
  JenkinsTraubSolver<Real>::init_next_h_poly()
  {
    const auto eps = Real{100} * s_eps;
    // Synthetic division of H by the quadratic 1, u, v
    NormalizationType type = none;
    this->remquo_quadratic(this->m_order - 1, this->m_u, this->m_v,
			   this->m_H, this->m_H_quot, this->m_c, this->m_d);
    if (std::abs(this->m_c) > eps * std::abs(this->m_H[this->m_order - 1])
     || std::abs(this->m_d) > eps * std::abs(this->m_H[this->m_order - 2]))
      {
	if (std::abs(this->m_d) < std::abs(this->m_c))
	  {
	    type = divide_by_c;
	    this->m_e = this->m_a / this->m_c;
	    this->m_f = this->m_d / this->m_c;
	    this->m_g = this->m_u * this->m_e;
	    this->m_h = this->m_v * this->m_b;
	    this->m_a3 = this->m_a * this->m_e
			+ (this->m_h / this->m_c + this->m_g) * this->m_b;
	    this->m_a1 = this->m_b - this->m_a * (this->m_d / this->m_c);
	    this->m_a7 = this->m_a
			+ this->m_g * this->m_d
			+ this->m_h * this->m_f;
	    return type;
	  }
	else
	  {
	    type = divide_by_d;
	    this->m_e = this->m_a / this->m_d;
	    this->m_f = this->m_c / this->m_d;
	    this->m_g = this->m_u * this->m_b;
	    this->m_h = this->m_v * this->m_b;
	    this->m_a3 = (this->m_a + this->m_g) * this->m_e
			+ this->m_h * (this->m_b / this->m_d);
	    this->m_a1 = this->m_b * this->m_f - this->m_a;
	    this->m_a7 = (this->m_f + this->m_u) * this->m_a + this->m_h;
	    return type;
	  }
      }
    else
      {
	type = near_h_root;
	return type;
      }
  }


/**
 * Computes the next H polynomials using scalars computed in init_next_h_poly.
 */
template<typename Real>
  void
  JenkinsTraubSolver<Real>::next_h_poly(NormalizationType type)
  {
    if (type == near_h_root)
      {
	// Use unscaled form of the recurrence if type is 3.
	this->m_H[0] = Real{0};
	this->m_H[1] = Real{0};
	for (int i = 2; i < this->m_order; ++i)
	  this->m_H[i] = this->m_H_quot[i-2];
	return;
      }
    auto ab_temp = this->m_a;
    if (type == divide_by_c)
      ab_temp = this->m_b;
    if (std::abs(this->m_a1) <= std::abs(ab_temp) * s_eps * Real{10})
      {
	// If a1 is nearly zero then use a special form of the recurrence.
	this->m_H[0] = Real{0};
	this->m_H[1] = -this->m_a7 * this->m_P_quot[0];
	for(int i = 2; i < this->m_order; ++i)
	  this->m_H[i] = this->m_a3 * this->m_H_quot[i-2]
			- this->m_a7 * this->m_P_quot[i-1];
      }
    else
      {
	// Use scaled form of the recurrence.
	this->m_a7 /= this->m_a1;
	this->m_a3 /= this->m_a1;
	this->m_H[0] = this->m_P_quot[0];
	this->m_H[1] = this->m_P_quot[1] - this->m_a7 * this->m_P_quot[0];
	for (int i = 2; i < this->m_order; ++i)
	  this->m_H[i] = this->m_a3 * this->m_H_quot[i-2]
			- this->m_a7 * this->m_P_quot[i-1]
			+ this->m_P_quot[i];
      }
  }

/**
 * Compute new estimates of the quadratic coefficients
 * using the scalars computed in init_next_h_poly.
 */
template<typename Real>
  std::pair<Real, Real>
  JenkinsTraubSolver<Real>::quadratic_coefficients(NormalizationType type)
  {
    if (type == near_h_root)
      return std::make_pair(Real{0}, Real{0});

    Real a4, a5;
    if (type == divide_by_d)
      {
	a4 = (this->m_a + this->m_g) * this->m_f + this->m_h;
	a5 = (this->m_f + this->m_u) * this->m_c + this->m_v * this->m_d;
      }
    else
      {
	a4 = this->m_a + this->m_u * this->m_b + this->m_h * this->m_f;
	a5 = this->m_c + (this->m_u + this->m_v * this->m_f) * this->m_d;
      }

    // Evaluate new quadratic coefficients.
    const auto n = this->m_order;
    const auto b1 = -this->m_H[n - 1] / this->m_P[n];
    const auto b2 = -(this->m_H[n - 2] + b1 * this->m_P[n - 1]) / this->m_P[n];
    const auto c1 = this->m_v * b2 * this->m_a1;
    const auto c2 = b1 * this->m_a7;
    const auto c3 = b1 * b1 * this->m_a3;
    const auto c4 = c1 - c2 - c3;
    if (auto temp = a5 + b1 * a4 - c4; temp == Real{0})
      return std::make_pair(Real{0}, Real{0});
    else
      {
	auto uu = this->m_u - (this->m_u * (c3 + c2)
		+ this->m_v * (b1 * this->m_a1 + b2 * this->m_a7)) / temp;
	auto vv = this->m_v * (Real{1} + c4 / temp);
	return std::make_pair(uu, vv);
      }
  }

/**
 * Divides the polynomial P by the quadratic 1x^2 + ux + v
 * placing the quotient in q and the remainder in a, b.
 */
template<typename Real>
  void
  JenkinsTraubSolver<Real>::remquo_quadratic(int nn, Real u, Real v,
					       std::vector<Real>& poly,
					       std::vector<Real>& quot,
					       Real& a, Real& b)
  {
    b = poly[0];
    quot[0] = b;
    a = poly[1] - b * u;
    quot[1] = a;
    for (int i = 2; i <= nn; ++i)
      {
	auto c = poly[i] - a * u - b * v;
	quot[i] = c;
	b = a;
	a = c;
      }	
  }


/**
 * Calculate the zeros of the quadratic az^2 + bz + c.
 * The quadratic formula, modified to avoid overflow, is used
 * to find the larger zero if the zeros are real and both
 * are complex. The smaller real zero is found directly from
 * the product of the zeros c/a.
 */
template<typename Real>
  void
  JenkinsTraubSolver<Real>::quadratic(Real a, Real b, Real c,
					Solution<Real>& z_small,
					Solution<Real>& z_large)
  {
    z_small = {};
    z_large = {};
    if (a == Real{0})
      { // Less than two roots.
	if (b != Real{0})
	  z_small = -c / b;
	return;
      }

    if (c == Real{0})
      { // one real root, one zero root.
	z_small = Real{0};
	z_large = -b / a;
	return;
      }

    // Compute discriminant avoiding overflow.
    auto b2 = b / Real{2};

    Real d, e;
    if (std::abs(b2) < std::abs(c))
      {
	e = std::copysign(a, c);
	e = b2 * (b2 / std::abs(c)) - e;
	d = std::sqrt(std::abs(e)) * std::sqrt(std::abs(c));
      }
    else
      {
	e = Real{1} - (a / b2) * (c / b2);
	d = std::sqrt(std::abs(e)) * std::abs(b2);
      }

    if (e < Real{0})
      { // complex conjugate zeros.
	z_small = std::complex<Real>(-b2 / a, +std::abs(d / a));
	z_large = std::complex<Real>(-b2 / a, -std::abs(d / a));
      }
    else
      {
	if (b2 >= Real{0})
	  d = -d; // Real zeros.
	z_large = (-b2 + d) / a;
	z_small = Real{0};
	if (z_large != Real{0})
	  z_small = (c / z_large) / a;
      }
  }

} // namespace emsr

#endif // SOLVER_JENKINS_TRAUB_TCC
