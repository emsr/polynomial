/*
$HOME/bin/bin/g++ -std=c++2a -Iinclude -g -Wall -Wextra -o test_quadratic_step test_quadratic_step.cpp
*/

#include <limits>
#include <complex>
#include <ext/polynomial.h>

#include <iostream>

  /**
   * I think this is trying to factor out a quadratic
   * from a complex-coefficient polynomial.
   */
  template<typename _Tp>
    void
    __root_quadratic(__gnu_cxx::_Polynomial<std::complex<_Tp>>& __p,
		     std::complex<_Tp>& __b, std::complex<_Tp>& __c, _Tp __eps)
    {
      using _Poly = __gnu_cxx::_Polynomial<std::complex<_Tp>>;
      constexpr int _S_max_iter = 50;
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_tiny = _Tp{100} * _S_eps;
      //auto __n = __p.degree();
      _Poly __q, __qq, __rem;
      for (int __iter = 0; __iter < _S_max_iter; ++__iter)
	{
	  _Poly __d({__c, __b, _Tp{1}});

	  // First division: r, s.
	  divmod(__p, __d, __q, __rem);
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
	  auto __delc = (-__r * __sb + __s * __rb) * __dv;
	  __b += __delb;
	  __delc = (-__r * __sb + __s * __rb) * __dv;
	  __c += __delc;
	  if ((std::abs(__delb) <= __eps * std::abs(__b)
	      || std::abs(__b) < _S_tiny)
           && (std::abs(__delc) <= __eps * std::abs(__c)
	      || std::abs(__c) < _S_tiny))
	    return;
	}
      std::__throw_runtime_error(__N("qroot: "
				     "Maximum number of iterations exceeded"));
    }

template<typename _Tp>
  void
  test_quadratic_step()
  {
    constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();

    __gnu_cxx::_Polynomial<std::complex<_Tp>>
    CP({std::complex<_Tp>(0.0, -1.0),
	std::complex<_Tp>(1.0, -2.0),
	std::complex<_Tp>(2.0, -3.0),
	std::complex<_Tp>(3.0, -4.0)});

    std::cout << "\ntest_quadratic\n";
    std::cout << "CP = " << CP << '\n';

    for (unsigned i = 0; i <= CP.degree(); i += 2)
      {
	std::cout << '\n';
	std::complex<_Tp> b, c;
	const auto eps = 100 * std::numeric_limits<_Tp>::epsilon();
	__root_quadratic(CP, b, c, eps);
	std::cout << "b   = " << b << '\n';
	std::cout << "c   = " << c << '\n';
	//std::cout << "CP(x) = " << CP(x) << '\n';
	__gnu_cxx::_Polynomial<std::complex<_Tp>>
	  CZ({c, b, std::complex<_Tp>(1)});
	CP.deflate(CZ, 10 * _S_eps);
	std::cout << "CP = " << CP << '\n';
      }
  }

int
main()
{
  test_quadratic_step<double>();
}
