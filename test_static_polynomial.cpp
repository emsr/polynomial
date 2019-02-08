/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_static_polynomial test_static_polynomial.cpp -lquadmath
./test_static_polynomial

g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -o test_static_polynomial test_static_polynomial.cpp -lquadmath
*/

//  Get past a bug....
// $HOME/bin/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__=0 -o test_static_polynomial test_static_polynomial.cpp

#include <initializer_list>
#include <iostream>
#include <complex>
#include <sstream>

#include "ext/static_polynomial.h"

constexpr void
test_static_polynomial()
{
}

int
main()
{
  float aa[5]{1.0f, 2.0f, -1.5f, 0.2f, -0.1f};

  __gnu_cxx::_StaticPolynomial<float, 5> p(aa);
  std::cout << "p: " << p << '\n';
  auto yp = p(1.2);
  std::cout << "p(1.2) = " << yp << '\n';

  __gnu_cxx::_StaticPolynomial<float, 5> a{{1.0f, 2.0f, -1.5f, 0.2f, -0.1f}};
  std::cout << "a: " << a << '\n';
  auto ya = a(1.2);
  std::cout << "a(1.2) = " << ya << '\n';

  auto za = __gnu_cxx::_StaticPolynomial<float, 5>({1.0f, 2.0f, -1.5f, 0.2f, -0.1f})(1.2);
  std::cout << "z(1.2) = " << za << '\n';

  double cc[1]{1.0};
  __gnu_cxx::_StaticPolynomial<double, 1> c(cc);
  std::cout << "c: " << c << '\n';

  auto dc = c.derivative();
  std::cout << "dc: " << dc << '\n';

  auto ic = c.integral(3.0);
  std::cout << "ic: " << ic << '\n';
}
