#ifndef SOLVER_QUADRATIC_H
#define SOLVER_QUADRATIC_H 1


#include <complex>

#include <ext/polynomial.h>
//#include <ext/solution.h> // For solution_t

namespace __gnu_cxx
{

  /**
   * A solver for complex-coefficient polynomials due to Laguerre.
   */
  template<typename _Real>
    class _QuadraticSolver
    {
    public:

      _QuadraticSolver(_Polynomial<std::complex<_Real>>& _P)
      : _M_poly(_P), _M_num_iters{0}
      { }

      //std::vector<solution_t<_Real>> solve();
      std::vector<std::complex<_Real>> solve();

      _Polynomial<std::complex<_Real>>
      step()
      {
	const auto __q = this->_M_root_quadratic();
	this->_M_poly.deflate(__q, _Real{10} * _S_eps);
	return __q;
      }

      int
      num_iters() const
      { return this->_M_num_iters; }

      int
      max_num_iters() const
      { return this->_M_max_num_iters; }

      const _Polynomial<std::complex<_Real>>&
      polynomial() const
      { return this->_M_poly; }

    private:

      // Estimated fractional roundoff error.
      static constexpr _Real _S_eps = std::numeric_limits<_Real>::epsilon();
      static constexpr _Real _S_tiny = _Real{10} * _S_eps;

      // Fractional roundoff error.
      _Real _M_eps = _Real{100} * _S_eps;

      int _M_max_iter = 50;

      _Polynomial<std::complex<_Real>> _M_root_quadratic();

      _Polynomial<std::complex<_Real>> _M_poly;

      int _M_num_iters = 0;
    };


} // namespace __gnu_cxx

#include <ext/solver_quadratic.tcc>

#endif // SOLVER_QUADRATIC_H
