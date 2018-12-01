/*
$HOME/bin/bin/g++ -std=c++2a -Iinclude -g -Wall -Wextra -o test_laguerre_step test_laguerre_step.cpp
*/

#include <limits>
#include <complex>
#include <ext/polynomial.h>

#include <iostream>

  /**
   * Find a root of a complex-coefficient polynomial by Laguerre's method.
   * This routine can be iterated by dividing the original polynomial
   * by the root factor (z - x) where x is the found root and finding
   * the next root.
   */
  template<typename _Tp>
    void
    __root_laguerre(__gnu_cxx::_Polynomial<std::complex<_Tp>>& __a,
		    std::complex<_Tp>& __x, int& __its)
    {
      // Estimated fractional roundoff error.
      constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();

      // Number of fractional values.
      constexpr int MR = 8;
      // Fractions used to break a limit cycle.
      static const _Tp
      _S_frac[MR + 1]
      {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};

      // Number of steps taken before trying a new fraction.
      constexpr int MT = 10;

      constexpr int _S_max_iter = MT * MR;

      int __m = __a.degree();
      for (int __iter = 1; __iter <= _S_max_iter; ++__iter)
	{
	  __its = __iter;
	  // Efficient computation of the polynomial
	  // and its first two derivatives. f stores P''(x)/2.
	  auto __b = __a[__m];
	  auto __err = std::abs(__b);
	  const auto __abx = std::abs(__x);
	  std::complex<_Tp> __d{}, __f{};
	  for (int __j = __m - 1; __j >= 0; --__j)
	    {
	      __f = __x * __f + __d;
	      __d = __x * __d + __b;
	      __b = __x * __b + __a[__j];
	      __err = __abx * __err + std::abs(__b);
	    }
	  __err *= _S_eps;
	  // Estimate of roundoff error in evaluating polynomial.
	  if (std::abs(__b) <= __err) // We have the root.
	    return;
	  // Use Laguerre's formula.
	  const auto __g = __d / __b;
	  const auto __g2 = __g * __g;
	  const auto __h = __g2 - _Tp{2} * __f / __b;
	  const auto __sq = std::sqrt(_Tp(__m - 1) * (_Tp(__m) * __h - __g2));
	  auto __gp = __g + __sq;
	  const auto __gm = __g - __sq;
	  const auto __abp = std::abs(__gp);
	  const auto __abm = std::abs(__gm);
	  if (__abp < __abm)
	    __gp = __gm;
	  const auto __dx = std::max(__abp, __abm) > _Tp{0}
			  ? _Tp(__m) / __gp
			  : std::polar(_Tp{1} + __abx, _Tp(__iter));
	  const auto __x1 = __x - __dx;
	  if (__x == __x1)
	    return;
	  if (__iter % MT != 0)
	    __x = __x1;
	  else
	    __x -= _S_frac[__iter / MT] * __dx;
	}
      std::__throw_runtime_error(__N("__root_laguerre: "
				     "Maximum number of iterations exceeded"));
    }

template<typename _Tp>
  void
  test_laguerre_step()
  {
    __gnu_cxx::_Polynomial<std::complex<_Tp>>
    CP({std::complex<_Tp>(0.0, -1.0),
	std::complex<_Tp>(1.0, -2.0),
	std::complex<_Tp>(2.0, -3.0),
	std::complex<_Tp>(3.0, -4.0)});

    std::cout << "\ntest_laguerre\n";
    std::cout << "CP = " << CP << '\n';

    for (unsigned i = 0; i <= CP.degree(); ++i)
      {
	std::cout << '\n';
	std::complex<_Tp> x;
	int its = 0;
	__root_laguerre(CP, x, its);
	std::cout << "x   = " << x << '\n';
	std::cout << "its = " << its << '\n';
	std::cout << "CP(x) = " << CP(x) << '\n';
	__gnu_cxx::_Polynomial<std::complex<_Tp>>
	  CZ({x, std::complex<_Tp>(1)});
	CP /= CZ;
	std::cout << "CP = " << CP << '\n';
      }
  }

int
main()
{
  test_laguerre_step<double>();
}
