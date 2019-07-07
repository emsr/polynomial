// Math extensions -*- C++ -*-

// Copyright (C) 2018-2019 Free Software Foundation, Inc.
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

#include <ext/solver_low_degree.h>

namespace __gnu_cxx
{

/**
 * Constructor from input polynomial.
 */
template<typename _Real>
  _JenkinsTraubSolver<_Real>::
  _JenkinsTraubSolver(const std::vector<_Real>& __op)
  : _P(__op)
  {
    if (this->_P.size() == 0)
      std::__throw_domain_error("Polynomial degree must be at least 1.");

    // Algorithm fails of the leading coefficient is zero.
    // We could erase leading-order zero coefficients.
    if (this->_P[0] == _Real{0})
      std::__throw_domain_error("Leading coefficient must be nonzero.");

    const auto __degree = this->_P.size() - 1;
    this->__order = __degree;
    this->_P_quot.resize(__degree + 1);
    this->_H.resize(__degree + 1);
    this->_H_quot.resize(__degree + 1);
    this->_H_save.resize(__degree + 1);
  }

/**
 *
 */
template<typename _Real>
  std::vector<solution_t<_Real>>
  _JenkinsTraubSolver<_Real>::solve()
  {
    // Initialization of constants for shift rotation.
    auto __xx = __gnu_cxx::math::__one_div_root_2_v<_Real>;
    auto __yy = -__xx;
    const auto __cosr = std::cos(_S_rotation);
    const auto __sinr = std::sin(_S_rotation);

    std::vector<solution_t<_Real>> __zero;
    __zero.reserve(this->_P.size());

    // Remove the zeros at the origin, if any.
    while (this->_P[this->__order] == _Real{0})
      {
	__zero.push_back(_Real{0});
	--this->__order;
      }
    if (this->__order < 1)
      return __zero;

    std::vector<_Real> __pt(this->__order + 1);
    std::vector<_Real> _H_temp(this->__order + 1);

    while (true)
      {
	// Start the algorithm for one zero.
	this->__num_iters = 0;
	if (this->__order == 1)
	  {
	    __zero.push_back(-this->_P[1] / this->_P[0]);
	    --this->__order;
	    return __zero;
	  }
	// Calculate the final zero or pair of zeros.
	if (this->__order == 2)
	  {
	    solution_t<_Real> __z_small, __z_large;
	    this->quadratic(this->_P[0], this->_P[1], this->_P[2],
			    __z_small, __z_large);
	    if (__z_small.index() != 0)
	      {
		__zero.push_back(__z_small);
		--this->__order;
	      }
	    if (__z_large.index() != 0)
	      {
		__zero.push_back(__z_large);
		--this->__order;
	      }
	    return __zero;
	  }

	// Find largest and smallest moduli of coefficients.
	auto __a_max = _Real{0};
	auto __a_min = _S_huge;
	for (int __i = 0; __i <= this->__order; ++__i)
	  {
	    const auto __x = std::abs(this->_P[__i]);
	    if (__x > __a_max)
	      __a_max = __x;
	    if (__x != _Real{0} && __x < __a_min)
	      __a_min = __x;
	  }
	// Scale if there are large or very tiny coefficients.
	// Computes a scale factor to multiply the coefficients
	// of the polynomial. The scaling is done to avoid overflow
	// and to avoid undetected underflow interfering
	// with the convergence criterion.
	// The factor is a power of the base.
	auto __scale = _S_low / __a_min;
	bool __rescale = true;
	if (__scale > _Real{1} && _S_huge / __scale < __a_max)
	  __rescale = false;
	if (__scale <= _Real{1} && __a_max < _Real{10})
	  __rescale = false;

	if (__rescale)
	  {
	    // Scale polynomial.
	    if (__scale == _Real{0})
	      __scale = _S_tiny;
	    const auto __l = std::ilogb(__scale);
	    const auto __factor = std::pow(_S_base, __l);
	    if (__factor != _Real{1})
	      for (int __i = 0; __i <= this->__order; ++__i)
		this->_P[__i] *= __factor;
	  }

	// Compute lower bound on moduli of roots.
	for (int __i = 0; __i <= this->__order; ++__i)
	  __pt[__i] = std::abs(this->_P[__i]);
	__pt[this->__order] = -__pt[this->__order];
	// Compute upper estimate of bound.
	auto __x = std::exp((std::log(-__pt[this->__order])
			    - std::log(__pt[0])) / _Real(this->__order));
	// If Newton step at the origin is better, use it.	
	if (__pt[this->__order - 1] != _Real{0})
	  {
	    const auto __xm = -__pt[this->__order] / __pt[this->__order - 1];
	    if (__xm < __x)
	      __x = __xm;
	  }
	// Chop the interval (0,x) until ff <= 0.
	while (true)
	  {
	    auto __xm = __x * _Real{0.1L};
	    auto __ff = __pt[0];
	    for (int __i = 1; __i <= this->__order; ++__i)
	      __ff = __ff * __xm + __pt[__i];
	    if (__ff <= _Real{0})
	      break;
	    __x = __xm;
	  }
	// Do Newton interation until x converges to two decimal places.
	auto __dx = __x;
	while (std::abs(__dx / __x) > this->min_log_deriv)
	  {
	    auto __ff = __pt[0];
	    auto __df = __ff;
	    for (int __i = 1; __i < this->__order; ++__i)
	      {
		__ff = __ff * __x + __pt[__i];
		__df = __df * __x + __ff;
	      }
	    __ff = __ff * __x + __pt[this->__order];
	    __dx = __ff / __df;
	    __x -= __dx;
	    ++this->__num_iters;
	  }
	const auto bound = __x;
	// Compute the derivative as the initial _H polynomial
	// and do 5 steps with no shift.
	const auto __nm1 = this->__order - 1;
	for (int __i = 1; __i < this->__order; ++__i)
	  this->_H[__i] = _Real(this->__order - __i) * this->_P[__i]
			/ _Real(this->__order);
	this->_H[0] = this->_P[0];
	const auto __aa = this->_P[this->__order];
	const auto __bb = this->_P[this->__order - 1];
	this->__zerok = (this->_H[this->__order - 1] == _Real{0});
	for(int __jj = 0; __jj < 5; ++__jj)
	  {
	    ++this->__num_iters;
	    auto __cc = this->_H[this->__order - 1];
	    if (!this->__zerok)
	      {
		// Use a scaled form of recurrence if value of H at 0
		// is nonzero.	
		const auto __t = -__aa / __cc;
		for (int __i = 0; __i < __nm1; ++__i)
		  {
		    const auto __j = this->__order - __i - 1;
		    this->_H[__j] = __t * this->_H[__j - 1] + this->_P[__j];
		  }
		this->_H[0] = this->_P[0];
		this->__zerok = (std::abs(this->_H[this->__order - 1])
			    <= _Real{10} * _S_eps * std::abs(__bb));
	    }
	    else
	      {
		// Use unscaled form of recurrence.
		for (int __i = 0; __i < __nm1; ++__i)
		  {
		    const auto __j = this->__order - __i - 1;
		    this->_H[__j] = this->_H[__j - 1];
		  }
		this->_H[0] = _Real{0};
		this->__zerok = (this->_H[this->__order - 1] == _Real{0});
	      }
	  }
	// Save H for restarts with new shifts.
	_H_temp = this->_H;

	// Loop to select the quadratic corresponding to each new shift.
	for (int __count = 0; __count < 20; ++__count)
	  {
	    /*  Quadratic corresponds to a _Real shift to a	
	     *  non-real point and its complex conjugate. The point
	     *  has modulus bound and amplitude rotated by 94 degrees
	     *  from the previous shift.
	     */
	    const auto __xxx = __cosr * __xx - __sinr * __yy;
	    __yy = __sinr * __xx + __cosr * __yy;
	    auto __xx = __xxx;
	    this->__sr = bound * __xx;
	    this->__si = bound * __yy;
	    this->__u = -_Real{2} * this->__sr;
	    this->__v = bound;
	    auto __num_zeros = this->fxshfr(20 * (__count + 1));
	    bool __cycle = false;
	    if (__num_zeros != 0)
	      {
	      /*  The second stage jumps directly to one of the third
	       *  stage iterations and returns here if successful.
	       *  Deflate the polynomial, store the zero or zeros
	       *  and return to the main algorithm.
	       */
		__zero.push_back(this->__z_small);
		this->__order -= __num_zeros;
		this->_P = this->_P_quot;
		if (__num_zeros != 1)
		  __zero.push_back(this->__z_large);
		__cycle = true;
		break;
	      }
	    if (__cycle)
	      continue;

	    // If the iteration is unsuccessful another quadratic
	    // is chosen after restoring H.
	    this->_H = _H_temp;
	 }
      }
  }


/**
 * Computes up to L2 fixed shift H-polynomials, testing for convergence
 * in the linear or quadratic case.
 * Initiates one of the variable shift iterations and returns
 * the number of zeros found.
 */
template<typename _Real>
  int
  _JenkinsTraubSolver<_Real>::fxshfr(int __l2)
  {
    _Real __ots, __otv;
    int __iflag;

    int __num_zeros = 0;

    auto __betav = _Real{0.25};
    auto __betas = _Real{0.25};
    auto __oss = this->__sr;
    auto __ovv = this->__v;
    // Evaluate polynomial by synthetic division.
    this->remquo_quadratic(this->__order, this->__u, this->__v,
			   this->_P, this->_P_quot,
			   this->__a, this->__b);
    auto __type = this->init_next_h_poly();
    for (int __j = 0; __j < __l2; ++__j)
      {
	// Calculate next H polynomial and estimate v.
	this->next_h_poly(__type);
	__type = this->init_next_h_poly();
	auto [__ui, __vi] = this->quadratic_coefficients(__type);
	auto __vv = __vi;
	// Estimate s.
	auto __ss = _Real{0};
	if (this->_H[this->__order - 1] != _Real{0})
	  __ss = -this->_P[this->__order] / this->_H[this->__order - 1];
	auto __tv = _Real{1};
	auto __ts = _Real{1};
	if (__j == 0 || __type == near_h_root)
	  {
	    __ovv = __vv;
	    __oss = __ss;
	    __otv = __tv;
	    __ots = __ts;
	    continue;
	  }

	// Compute relative measures of convergence of s and v sequences.
	if (__vv != _Real{0})
	  __tv = std::abs((__vv - __ovv) / __vv);
	if (__ss != _Real{0})
	  __ts = std::abs((__ss - __oss) / __ss);

	// If decreasing, multiply two most recent convergence measures.
	const auto __tvv = __tv < __otv ? __tv * __otv : _Real{1};
	const auto __tss = __ts < __ots ? __ts * __ots : _Real{1};

	// Compare with convergence criteria.
	const auto __vpass = __tvv < __betav;
	const auto __spass = __tss < __betas;
	if (!(__spass || __vpass))
	  {
	    __ovv = __vv;
	    __oss = __ss;
	    __otv = __tv;
	    __ots = __ts;
	    continue;
	  }

	// At least one sequence has passed the convergence test.
	// Store variables before iterating.
	const auto __u_save = this->__u;
	const auto __v_save = this->__v;
	this->_H_save = this->_H;
	const auto __s = __ss;

	// Choose iteration according to the fastest converging sequence.
	auto __vtry = false;
	auto __stry = false;
	if ((__spass && !__vpass) || __tss < __tvv)
	  goto _TRY_LINEAR_;

  _TRY_QUADRATIC_:
	__num_zeros = this->iter_quadratic(__ui, __vi);
	if (__num_zeros > 0)
	  return __num_zeros;
	// Quadratic iteration has failed. Flag that it has
	// been tried and decrease the convergence criterion.
	__vtry = true;
	__betav *= _Real{0.25};
	// Try linear iteration if it has not been tried and
	// the S sequence is converging.
	if (__stry || !__spass)
	  goto _RESTORE_VARS_;
	this->_H = this->_H_save;

  _TRY_LINEAR_:
	__num_zeros = this->iter_real(__s, __iflag);
	if (__num_zeros > 0)
	  return __num_zeros;
	// Linear iteration has failed. Flag that it has been
	// tried and decrease the convergence criterion.
	__stry = true;
	__betas *= _Real{0.25};
	if (__iflag == 0)
	  goto _RESTORE_VARS_;
	// If linear iteration signals an almost real
	// zero attempt quadratic iteration.
	__ui = -_Real{2} * __s;
	__vi = __s * __s;
	goto _TRY_QUADRATIC_;

  _RESTORE_VARS_:
	// Restore variables.
	this->__u = __u_save;
	this->__v = __v_save;
	this->_H = this->_H_save;

	// Try quadratic iteration if it has not been tried
	// and the V sequence is converging.
	if (!__vtry && __vpass)
	  goto _TRY_QUADRATIC_;

	// Recompute polynomial quotient and remainder
        // to continue the second stage.
	this->remquo_quadratic(this->__order, this->__u, this->__v,
			       this->_P, this->_P_quot,
			       this->__a, this->__b);
	__type = this->init_next_h_poly();
      }
    return __num_zeros;
  }


/**
 * Variable-shift H-polynomial iteration for a quadratic factor
 * converges only if the zeros are equimodular or nearly so.
 * @param uu The linear coefficient of the starting quadratic equation.
 * @param vv The constant  coefficient of the starting quadratic equation.
 * @return The number of zeros found.
 */
template<typename _Real>
  int
  _JenkinsTraubSolver<_Real>::iter_quadratic(_Real __uu, _Real __vv)
  {
    _Real __mp, __omp, __ee, __relstp, __t, __zm;
    NormalizationType __type;

    int __num_zeros = 0;
    bool __tried = false;
    this->__u = __uu;
    this->__v = __vv;
    int __j = 0;

    while (true)
      {
	++this->__num_iters;
	this->quadratic(_Real{1}, this->__u, this->__v,
			this->__z_small, this->__z_large);
	// Return if roots of the quadratic are real and not
	// close to multiple or nearly equal and of opposite sign.
	if (std::abs(std::abs(real(this->__z_small))
		   - std::abs(real(this->__z_large)))
	       > _Real{0.01L} * std::abs(real(this->__z_large)))
	  return __num_zeros;
	// Evaluate polynomial by quadratic synthetic division.
	this->remquo_quadratic(this->__order, this->__u, this->__v,
			       this->_P, this->_P_quot, this->__a, this->__b);
	__mp = std::abs(this->__a - real(this->__z_small) * this->__b)
	   + std::abs(imag(this->__z_small) * this->__b);
	// Compute a rigorous bound on the rounding error in evaluating _P.
	__zm = std::sqrt(std::abs(this->__v));
	__ee = _Real{2} * std::abs(this->_P_quot[0]);
	__t = -real(this->__z_small) * this->__b;
	for (int __i = 1; __i < this->__order; ++__i)
	  __ee = __ee * __zm + std::abs(this->_P_quot[__i]);
	__ee = __ee * __zm + std::abs(this->__a + __t);
	__ee *= (_Real{5} * this->__mre + _Real{4} * this->__are);
	__ee -= (_Real{5} * this->__mre + _Real{2} * this->__are)
	    * (std::abs(this->__a + __t) + std::abs(this->__b) * __zm);
	__ee += _Real{2} * this->__are * std::abs(__t);
	// Iteration has converged sufficiently if the
	// polynomial value is less than 20 times this bound.
	if (__mp <= _Real{20} * __ee)
	  {
	    __num_zeros = 2;
	    return __num_zeros;
	  }
	++__j;
	// Stop iteration after 20 steps.
	if (__j > this->max_iter_quadratic)
	  return __num_zeros;
	if (__j < 2 || __relstp > _Real{0.01L} || __mp < __omp || __tried)
	  {
	    __omp = __mp;
	    // Calculate next H polynomial and new u and v.
	    __type = this->init_next_h_poly();
	    this->next_h_poly(__type);
	    __type = this->init_next_h_poly();
	    const auto [__ui, __vi] = this->quadratic_coefficients(__type);
	    // If vi is zero the iteration is not converging.
	    if (__vi == _Real{0})
	      return __num_zeros;
	    __relstp = std::abs((__vi - this->__v) / __vi);
	    this->__u = __ui;
	    this->__v = __vi;
	    continue;
	  }
	// A cluster appears to be stalling the convergence.
	// Five fixed shift steps are taken with a u, v close to the cluster.
	if (__relstp < _S_eps)
	  __relstp = _S_eps;
	__relstp = std::sqrt(__relstp);
	this->__u -= this->__u * __relstp;
	this->__v += this->__v * __relstp;
	this->remquo_quadratic(this->__order, this->__u, this->__v,
			       this->_P, this->_P_quot,
			       this->__a, this->__b);
	for (int __i = 0; __i < 5; ++__i)
	  {
	    __type = this->init_next_h_poly();
	    this->next_h_poly(__type);
	  }
	__tried = true;
	__j = 0;
      }
  }


/**
 * Variable-shift H polynomial iteration for a real zero.
 * @param sss Starting iterate.
 * @param iflag Flag to indicate a pair of zeros near real axis.
 * @return The number of zeros found.
 */
template<typename _Real>
  int
  _JenkinsTraubSolver<_Real>::iter_real(_Real __sss, int& __iflag)
  {
    auto __t = _Real{0};
    decltype(std::abs(this->_P[0])) __omp;

    int __num_zeros = 0;
    auto __s = __sss;
    __iflag = 0;
    int i_real = 0;

    while (true)
      {
	++this->__num_iters;
	auto __pval = this->_P[0];
	// Evaluate P at s.
	this->_P_quot[0] = __pval;
	for (int __i = 1; __i <= this->__order; ++__i)
	  {
	    __pval = __pval * __s + this->_P[__i];
	    this->_P_quot[__i] = __pval;
	  }
	auto __mp = std::abs(__pval);
	// Compute a rigorous bound on the error in evaluating P.
	const auto __ms = std::abs(__s);
	auto __ee = (this->__mre / (this->__are + this->__mre))
		  * std::abs(this->_P_quot[0]);
	for (int __i = 1; __i <= this->__order; ++__i)
	  __ee = __ee * __ms + std::abs(this->_P_quot[__i]);
	// Iteration has converged sufficiently if the polynomial
	// value is less than 20 times this bound.
	if (__mp <= _Real{20}
		 * ((this->__are + this->__mre) * __ee - this->__mre * __mp))
	  {
	    __num_zeros = 1;
	    this->__z_small = __s;
	    return __num_zeros;
	  }
	++i_real;
	// Stop iteration after max_iter_real steps.
	if (i_real > this->max_iter_real)
	  return __num_zeros;
	else if (i_real < 2
	  || std::abs(__t) > _Real{0.001L} * std::abs(__s - __t)
	  || __mp < __omp)
	  {
	    // Return if the polynomial value has increased significantly.
	    __omp = __mp;

	    // Compute t, the next polynomial, and the new iterate.
	    auto __hval = this->_H[0];
	    this->_H_quot[0] = __hval;
	    for (int __i = 1; __i < this->__order; ++__i)
	      {
		__hval = __hval * __s + this->_H[__i];
		this->_H_quot[__i] = __hval;
	      }

	    if (std::abs(__hval)
		 <= std::abs(this->_H[this->__order - 1]) * _Real{10} * _S_eps)
	      {
		// Use unscaled form.
		this->_H[0] = _Real{0};
		for (int __i = 1; __i < this->__order; ++__i)
		  this->_H[__i] = this->_H_quot[__i-1];
	      }
	    else
	      {
		// Use the scaled form of the recurrence if the value
		// of H at s is nonzero.
		__t = -__pval / __hval;
		this->_H[0] = this->_P_quot[0];
		for (int __i = 1; __i < this->__order; ++__i)
		  this->_H[__i] = __t * this->_H_quot[__i - 1]
				+ this->_P_quot[__i];
	      }

	    __hval = this->_H[0];
	    for (int __i = 1; __i < this->__order; ++__i)
	      __hval = __hval * __s + this->_H[__i];
	    auto __t = _Real{0};
	    if (std::abs(__hval)
		 > std::abs(this->_H[this->__order - 1] * _Real{10} * _S_eps))
	      __t = -__pval / __hval;
	    __s += __t;
	  }
	else
	  {
	    // A cluster of zeros near the real axis has been encountered.
	    // Return with iflag set to initiate a quadratic iteration.
	    __iflag = 1;
	    __sss = __s;
	    return __num_zeros;
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
template<typename _Real>
  typename _JenkinsTraubSolver<_Real>::NormalizationType
  _JenkinsTraubSolver<_Real>::init_next_h_poly()
  {
    const auto eps = _Real{100} * _S_eps;
    // Synthetic division of H by the quadratic 1, u, v
    NormalizationType __type = none;
    this->remquo_quadratic(this->__order - 1, this->__u, this->__v,
			   this->_H, this->_H_quot, this->__c, this->__d);
    if (std::abs(this->__c) > eps * std::abs(this->_H[this->__order - 1])
     || std::abs(this->__d) > eps * std::abs(this->_H[this->__order - 2]))
      {
	if (std::abs(this->__d) < std::abs(this->__c))
	  {
	    __type = divide_by_c;
	    this->__e = this->__a / this->__c;
	    this->__f = this->__d / this->__c;
	    this->__g = this->__u * this->__e;
	    this->__h = this->__v * this->__b;
	    this->__a3 = this->__a * this->__e
			+ (this->__h / this->__c + this->__g) * this->__b;
	    this->__a1 = this->__b - this->__a * (this->__d / this->__c);
	    this->__a7 = this->__a
			+ this->__g * this->__d
			+ this->__h * this->__f;
	    return __type;
	  }
	else
	  {
	    __type = divide_by_d;
	    this->__e = this->__a / this->__d;
	    this->__f = this->__c / this->__d;
	    this->__g = this->__u * this->__b;
	    this->__h = this->__v * this->__b;
	    this->__a3 = (this->__a + this->__g) * this->__e
			+ this->__h * (this->__b / this->__d);
	    this->__a1 = this->__b * this->__f - this->__a;
	    this->__a7 = (this->__f + this->__u) * this->__a + this->__h;
	    return __type;
	  }
      }
    else
      {
	__type = near_h_root;
	return __type;
      }
  }


/**
 * Computes the next H polynomials using scalars computed in init_next_h_poly.
 */
template<typename _Real>
  void
  _JenkinsTraubSolver<_Real>::next_h_poly(NormalizationType __type)
  {
    if (__type == near_h_root)
      {
	// Use unscaled form of the recurrence if type is 3.
	this->_H[0] = _Real{0};
	this->_H[1] = _Real{0};
	for (int __i = 2; __i < this->__order; ++__i)
	  this->_H[__i] = this->_H_quot[__i-2];
	return;
      }
    auto __ab_temp = this->__a;
    if (__type == divide_by_c)
      __ab_temp = this->__b;
    if (std::abs(this->__a1) <= std::abs(__ab_temp) * _S_eps * _Real{10})
      {
	// If a1 is nearly zero then use a special form of the recurrence.
	this->_H[0] = _Real{0};
	this->_H[1] = -this->__a7 * this->_P_quot[0];
	for(int __i = 2; __i < this->__order; ++__i)
	  this->_H[__i] = this->__a3 * this->_H_quot[__i-2]
			- this->__a7 * this->_P_quot[__i-1];
      }
    else
      {
	// Use scaled form of the recurrence.
	this->__a7 /= this->__a1;
	this->__a3 /= this->__a1;
	this->_H[0] = this->_P_quot[0];
	this->_H[1] = this->_P_quot[1] - this->__a7 * this->_P_quot[0];
	for (int __i = 2; __i < this->__order; ++__i)
	  this->_H[__i] = this->__a3 * this->_H_quot[__i-2]
			- this->__a7 * this->_P_quot[__i-1]
			+ this->_P_quot[__i];
      }
  }

/**
 * Compute new estimates of the quadratic coefficients
 * using the scalars computed in init_next_h_poly.
 */
template<typename _Real>
  std::pair<_Real, _Real>
  _JenkinsTraubSolver<_Real>::quadratic_coefficients(NormalizationType __type)
  {
    if (__type == near_h_root)
      return std::make_pair(_Real{0}, _Real{0});

    _Real __a4, __a5;
    if (__type == divide_by_d)
      {
	__a4 = (this->__a + this->__g) * this->__f + this->__h;
	__a5 = (this->__f + this->__u) * this->__c + this->__v * this->__d;
      }
    else
      {
	__a4 = this->__a + this->__u * this->__b + this->__h * this->__f;
	__a5 = this->__c + (this->__u + this->__v * this->__f) * this->__d;
      }

    // Evaluate new quadratic coefficients.
    const auto n = this->__order;
    const auto __b1 = -this->_H[n - 1] / this->_P[n];
    const auto __b2 = -(this->_H[n - 2] + __b1 * this->_P[n - 1])
		    / this->_P[n];
    const auto __c1 = this->__v * __b2 * this->__a1;
    const auto __c2 = __b1 * this->__a7;
    const auto __c3 = __b1 * __b1 * this->__a3;
    const auto __c4 = __c1 - __c2 - __c3;
    if (auto __temp = __a5 + __b1 * __a4 - __c4; __temp == _Real{0})
      return std::make_pair(_Real{0}, _Real{0});
    else
      {
	auto __uu = this->__u - (this->__u * (__c3 + __c2)
		+ this->__v * (__b1 * this->__a1 + __b2 * this->__a7)) / __temp;
	auto __vv = this->__v * (_Real{1} + __c4 / __temp);
	return std::make_pair(__uu, __vv);
      }
  }

/**
 * Divides the polynomial P by the quadratic 1x^2 + ux + v
 * placing the quotient in q and the remainder in a, b.
 */
template<typename _Real>
  void
  _JenkinsTraubSolver<_Real>::remquo_quadratic(int __nn, _Real __u, _Real __v,
					       std::vector<_Real>& __poly,
					       std::vector<_Real>& __quot,
					       _Real& __a, _Real& __b)
  {
    __b = __poly[0];
    __quot[0] = __b;
    __a = __poly[1] - __b * __u;
    __quot[1] = __a;
    for (int __i = 2; __i <= __nn; ++__i)
      {
	auto __c = __poly[__i] - __a * __u - __b * __v;
	__quot[__i] = __c;
	__b = __a;
	__a = __c;
      }	
  }


/**
 * Calculate the zeros of the quadratic az^2 + bz + c.
 * The quadratic formula, modified to avoid overflow, is used
 * to find the larger zero if the zeros are real and both
 * are complex. The smaller real zero is found directly from
 * the product of the zeros c/a.
 */
template<typename _Real>
  void
  _JenkinsTraubSolver<_Real>::quadratic(_Real a, _Real b, _Real c,
					solution_t<_Real>& __z_small,
					solution_t<_Real>& __z_large)
  {
    __z_small = {};
    __z_large = {};
    if (a == _Real{0})
      { // Less than two roots.
	if (b != _Real{0})
	  __z_small = -c / b;
	return;
      }

    if (c == _Real{0})
      { // one real root, one zero root.
	__z_small = _Real{0};
	__z_large = -b / a;
	return;
      }

    // Compute discriminant avoiding overflow.
    auto __b2 = b / _Real{2};

    _Real d, e;
    if (std::abs(__b2) < std::abs(c))
      {
	e = std::copysign(a, c);
	e = __b2 * (__b2 / std::abs(c)) - e;
	d = std::sqrt(std::abs(e)) * std::sqrt(std::abs(c));
      }
    else
      {
	e = _Real{1} - (a / __b2) * (c / __b2);
	d = std::sqrt(std::abs(e)) * std::abs(__b2);
      }

    if (e < _Real{0})
      { // complex conjugate zeros.
	__z_small = std::complex<_Real>(-__b2 / a, +std::abs(d / a));
	__z_large = std::complex<_Real>(-__b2 / a, -std::abs(d / a));
      }
    else
      {
	if (__b2 >= _Real{0})
	  d = -d; // Real zeros.
	__z_large = (-__b2 + d) / a;
	__z_small = _Real{0};
	if (__z_large != _Real{0})
	  __z_small = (c / __z_large) / a;
      }
  }

} // namespace __gnu_cxx

#endif // SOLVER_JENKINS_TRAUB_TCC
