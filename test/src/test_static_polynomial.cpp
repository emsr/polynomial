
//  Get past a bug....
// $HOME/bin/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__=0 -o test_static_polynomial test_static_polynomial.cpp

#include <initializer_list>
#include <iostream>
#include <complex>
#include <sstream>

#include <ext/static_polynomial.h>

constexpr void
test_static_polynomial()
{
}

int
main()
{
  std::cout.setf(std::ios_base::boolalpha);
  std::cout.precision(std::numeric_limits<double>::digits10);

  constexpr float aa[5]{1.0f, 2.0f, -1.5f, 0.2f, -0.1f};

  constexpr __gnu_cxx::_StaticPolynomial<float, 5> p(aa);
  std::cout << "p: " << p << '\n';
  auto yp = p(1.2);
  std::cout << "p(1.2) = " << yp << '\n';

  constexpr __gnu_cxx::_StaticPolynomial<float, 5> q{{1.0f, 2.0f, -1.5f, 0.2f, -0.1f}};
  std::cout << "q: " << q << '\n';
  auto yq = q(1.2);
  std::cout << "q(1.2) = " << yq << '\n';

  constexpr auto za = __gnu_cxx::_StaticPolynomial<float, 5>({1.0f, 2.0f, -1.5f, 0.2f, -0.1f})(1.2);
  std::cout << "z(1.2) = " << za << '\n';

  constexpr double cc[1]{1.0};
  constexpr __gnu_cxx::_StaticPolynomial<double, 1> c(cc);
  std::cout << "c: " << c << '\n';

  auto dc = c.derivative();
  std::cout << "dc: " << dc << '\n';

  auto ic = c.integral(3.0);
  std::cout << "ic: " << ic << '\n';

  constexpr __gnu_cxx::_StaticPolynomial<double, 4> P({0.0, 1.0, 2.0, 3.0});
  std::cout << "P = " << P << '\n';
  std::cout << "+P = " << +P << '\n';
  std::cout << "-P = " << -P << '\n';
  std::cout << "P = " << P << '\n';
  std::cout << "degree(P) = " << P.degree() << '\n';

  constexpr __gnu_cxx::_StaticPolynomial<double, 2> Q({2.0, 1.0});
  std::cout << "Q = " << Q << '\n';
  std::cout << "degree(Q) = " << Q.degree() << '\n';
  std::cout << "P + Q = " << P + Q << '\n';
  std::cout << "P - Q = " << P - Q << '\n';
  //std::cout << "P * Q = " << P * Q << '\n';
  std::cout << "P / Q = " << P / Q << '\n';
  std::cout << "P % Q = " << P % Q << '\n';

  constexpr double b = 5.0;
  std::cout << "b = " << b << '\n';
  std::cout << "P + b = " << P + b << '\n';
  std::cout << "P - b = " << P - b << '\n';
  std::cout << "P * b = " << P * b << '\n';
  std::cout << "P / b = " << P / b << '\n';

  constexpr double a = 2.0;
  std::cout << "a = " << a << '\n';
  std::cout << "a + P = " << a + P << '\n';
  std::cout << "a - P = " << a - P << '\n';
  std::cout << "a * P = " << a * P << '\n';
}
