#ifndef SOLVER_LAGUERRE_TCC
#define SOLVER_LAGUERRE_TCC 1

namespace __gnu_cxx
{

  /**
   * Find a root of a complex-coefficient polynomial by Laguerre's method.
   * This routine can be iterated by dividing the original polynomial
   * by the root factor (z - x) where x is the found root and finding
   * the next root.
   */
  template<typename _Real>
    std::complex<_Real>
    _LaguerreSolver<_Real>::_M_root_laguerre()
    {
      using __cmplx = std::complex<_Real>;

      __cmplx __x;
      int __m = this->_M_poly.degree();

      //if (__m == 1)
	//return -this->_M_poly.cefficient(0) / this->_M_poly.cefficient(1);

      const int __max_iter = this->_M_max_iter();
      for (int __iter = 1; __iter <= __max_iter; ++__iter)
	{
	  this->_M_num_iters = __iter;
	  // Efficient computation of the polynomial
	  // and its first two derivatives. F stores P''(x)/2.
	  auto __b = this->_M_poly[__m];
	  auto __err = std::abs(__b);
	  const auto __abx = std::abs(__x);
	  __cmplx __d{}, __f{};
	  for (int __j = __m - 1; __j >= 0; --__j)
	    {
	      __f = __x * __f + __d;
	      __d = __x * __d + __b;
	      __b = __x * __b + this->_M_poly[__j];
	      __err = __abx * __err + std::abs(__b);
	    }
	  __err *= _S_eps;
	  // Estimate of roundoff error in evaluating polynomial.
	  if (std::abs(__b) <= __err) // We have the root.
	    return __x;

	  // Use Laguerre's formula.
	  const auto __g = __d / __b;
	  const auto __g2 = __g * __g;
	  const auto __h = __g2 - _Real{2} * __f / __b;
	  const auto __sq = std::sqrt(_Real(__m - 1)
				   * (_Real(__m) * __h - __g2));
	  auto __gp = __g + __sq;
	  const auto __gm = __g - __sq;
	  const auto __abp = std::abs(__gp);
	  const auto __abm = std::abs(__gm);
	  if (__abp < __abm)
	    __gp = __gm;
	  const auto __dx = std::max(__abp, __abm) > _Real{0}
			  ? _Real(__m) / __gp
			  : std::polar(_Real{1} + __abx, _Real(__iter));
	  const auto __x1 = __x - __dx;
	  if (__x == __x1)
	    return __x;
	  if (__iter % this->_M_steps_per_frac != 0)
	    __x = __x1;
	  else
	    __x -= _S_frac[__iter / this->_M_steps_per_frac] * __dx;
	}

      std::__throw_runtime_error(__N("_M_root_laguerre: "
				     "Maximum number of iterations exceeded"));
    }

  template<typename _Real>
    std::vector<std::complex<_Real>>
    _LaguerreSolver<_Real>::solve()
    {
      using __cmplx = std::complex<_Real>;

      std::vector<__cmplx> __roots;
      __roots.reserve(this->_M_poly.degree());
      for (unsigned __i = 0; __i < this->_M_poly.degree(); ++__i)
	{
	  const auto __z = this->_M_root_laguerre();

	  _Polynomial<__cmplx> __zpoly({__z, __cmplx{1}});
	  this->_M_poly /= __zpoly;

	  __roots.push_back(__z);
	}
      return __roots;
    }

} // namespace __gnu_cxx

#endif // SOLVER_LAGUERRE_TCC
