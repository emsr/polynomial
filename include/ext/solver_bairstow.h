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
