// Math extensions -*- C++ -*-

// Copyright (C) 2019 Free Software Foundation, Inc.
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
 * @file solver_quadratic.tcc Class declaration for the quadratic solver.
 */

/**
 * @def  SOLVER_QUADRATIC_TCC
 *
 * @brief  A guard for the _QuadraticSolver class header.
 */
#ifndef SOLVER_QUADRATIC_TCC
#define SOLVER_QUADRATIC_TCC 1

namespace __gnu_cxx
{

  /**
   * I think this is trying to factor out a quadratic
   * from a complex-coefficient polynomial.
   */
  template<typename _Tp>
    _Polynomial<std::complex<_Tp>>
    _QuadraticSolver<_Tp>::_M_root_quadratic()
    {
      using _Cmplx = std::complex<_Tp>;
      using _Poly = _Polynomial<_Cmplx>;

      if (this->_M_poly.degree() <= 2)
	return this->_M_poly;

      this->_M_num_iters = 0;

      _Cmplx __c, __b;
      _Poly __q, __qq, __rem;
      for (int __iter = 0; __iter < this->_M_max_iter; ++__iter)
	{
	  ++this->_M_num_iters;

	  _Poly __d({__c, __b, _Cmplx{1}});

	  // First division: r, s.
	  divmod(this->_M_poly, __d, __q, __rem);
	  const auto __s = __rem[0];
	  const auto __r = __rem[1];

	  // Second division: partial r, s with respect to c.
	  divmod(__q, __d, __qq, __rem);
	  const auto __sc = -__rem[0];
	  const auto __rc = -__rem[1];
	  const auto __sb = -__c * __rc;
	  const auto __rb = -__b * __rc + __sc;

	  // Solve 2x2 equation.
	  const auto __dv = _Tp{1} / (__sb * __rc - __sc * __rb);
	  const auto __delb = ( __r * __sc - __s * __rc) * __dv;
	  __b += __delb;
	  const auto __delc = (-__r * __sb + __s * __rb) * __dv;
	  __c += __delc;
	  if ((std::abs(__delb) <= this->_M_eps * std::abs(__b)
	      || std::abs(__b) < _S_tiny)
           && (std::abs(__delc) <= this->_M_eps * std::abs(__c)
	      || std::abs(__c) < _S_tiny))
	    return _Poly({__c, __b, _Cmplx{1}});
	}
      std::__throw_runtime_error(__N("_M_root_quadratic: "
				     "Maximum number of iterations exceeded"));
    }

} // namespace __gnu_cxx

#endif // SOLVER_QUADRATIC_TCC
