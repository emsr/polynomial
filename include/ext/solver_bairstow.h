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
 * @file solver_bairstow.h Class declaration for the Bairstow solver.
 */

/**
 * @def  SOLVER_BAIRSTOW_H
 *
 * @brief  A guard for the _BairstowSolver class header.
 */
#ifndef SOLVER_BAIRSTOW_H
#define SOLVER_BAIRSTOW_H 1

#include <ext/solver_low_degree.h>

namespace __gnu_cxx
{

///
/// @todo: If you *don't* reverse the input array you solve for 1/z.
/// @todo Take a max_error.  If this->_M_eps grows larger than max_error throw.
///
template<typename _Real>
  class _BairstowSolver
  {
  public:

    _BairstowSolver(const std::vector<_Real>& __coeff,
		    unsigned int __seed = std::random_device()())
    : _M_coeff{__coeff.rbegin(), __coeff.rend()},
      _M_b(__coeff.size()), _M_c(__coeff.size()),
      _M_order(__coeff.size() - 1),
      _M_urng(__seed), _M_pdf(_Real{0}, _Real{2})
    {
      if (this->_M_coeff.size() == 0)
	std::__throw_domain_error("_BairstowSolver: "
				  "Coefficient size must be nonzero.");

      if (this->_M_coeff[0] == _Real{0})
	std::__throw_domain_error("_BairstowSolver: "
				  "Leading-order coefficient must be nonzero.");

      const auto __scale = this->_M_coeff[0];
      for (int __i = 0; __i <= this->_M_order; ++__i)
	this->_M_coeff[__i] /= __scale;

      this->_M_zero.reserve(this->_M_coeff.size());
    }

    std::vector<solution_t<_Real>> solve();
    std::vector<_Real> equations() const;

  private:

    void _M_iterate();

    template<int _Index>
    static void _S_refine_quadratic_eqn(_Real& __dr, _Real& __r,
					_Real& __ds, _Real& __s,
					std::array<solution_t<_Real>, 2>& __w);

    void
    _M_add_zero(solution_t<_Real> __z)
    {
      this->_M_zero.push_back(__z);
      --this->_M_order;
    }

    static constexpr _Real _S_ratio
	= _Real{std::numeric_limits<_Real>::digits}
	/ std::numeric_limits<double>::digits;
    static constexpr int _S_max_rand_iter = 200 * _S_ratio;
    static constexpr int _S_max_error_iter = 500 * _S_ratio;
    static constexpr auto _S_eps_factor = _Real{10} * _S_ratio;
    static constexpr auto _S_eps
      = _S_eps_factor * std::numeric_limits<_Real>::epsilon();

    std::vector<_Real> _M_coeff;
    std::vector<_Real> _M_b;
    std::vector<_Real> _M_c;
    std::vector<solution_t<_Real>> _M_zero;
    _Real _M_eps = _S_eps;
    int _M_order;
    bool _M_precision_error = false;
    std::mt19937 _M_urng;
    std::uniform_real_distribution<_Real> _M_pdf;
  };

} // namespace __gnu_cxx

#include <ext/solver_bairstow.tcc>

#endif // SOLVER_BAIRSTOW_H
