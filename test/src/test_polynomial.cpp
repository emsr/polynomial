
#include <cmath>
#include <iostream>
#include <complex>
#include <sstream>

#include <ext/polynomial.h>


// Quick factorial impl.
template<typename Real>
  constexpr Real
  factorial(unsigned int __n)
  {
    if (__n == 0)
      return Real{1};
    else
      {
	auto __fact = Real{1};
	for (unsigned int __i = 1; __i <= __n; ++__i)
	  __fact *= Real(__i);
	return __fact;
      }
  }

int
main()
{
  using namespace std::literals::complex_literals;

  using Real = double;
  using Poly = emsr::Polynomial<Real>;
  using Cmplx = std::complex<Real>;
  using CPoly = emsr::Polynomial<Cmplx>;

  std::cout.setf(std::ios_base::boolalpha);
  std::cout.precision(std::numeric_limits<Real>::digits10);

  int num_errors = 0;

  // Polynomial properties...
  Poly P({Real{0}, Real{1}, Real{2}, Real{3}});
  std::cout << "\nP = " << P << '\n';
  std::cout << "+P = " << +P << '\n';
  std::cout << "-P = " << -P << '\n';
  std::cout << "P = " << P << '\n';
  std::cout << "degree(P) = " << P.degree() << '\n';
  num_errors += P.degree() != 3;

  // Polynomial maths...
  Poly Q({Real{2}, Real{1}});
  std::cout << "\nQ = " << Q << '\n';
  std::cout << "degree(Q) = " << Q.degree() << '\n';

  std::cout << "P + Q = " << P + Q << '\n';
  num_errors += (P + Q) != Poly({Real{2},Real{2},Real{2},Real{3}});

  std::cout << "P - Q = " << P - Q << '\n';
  num_errors += (P - Q) != Poly({Real{-2},Real{0},Real{2},Real{3}});

  std::cout << "P * Q = " << P * Q << '\n';
  num_errors += (P * Q) != Poly({Real{0},Real{2},Real{5},Real{8},Real{3}});

  std::cout << "P / Q = " << P / Q << '\n';
  num_errors += (P / Q) != Poly({Real{9},Real{-4},Real{3}});

  std::cout << "P % Q = " << P % Q << '\n';
  num_errors += (P % Q) != Poly({Real{-18}});

  // Post scalar...
  Real b = 5;
  std::cout << "\nb = " << b << '\n';

  std::cout << "P + b = " << P + b << '\n';
  num_errors += (P + b) != Poly({Real{5},Real{1},Real{2},Real{3}});

  std::cout << "P - b = " << P - b << '\n';
  num_errors += (P - b) != Poly({Real{-5},Real{1},Real{2},Real{3}});

  std::cout << "P * b = " << P * b << '\n';
  num_errors += (P * b) != Poly({Real{0},Real{5},Real{10},Real{15}});

  std::cout << "P / b = " << P / b << '\n';
  num_errors += (P / b) != Poly({Real{0},Real{0.2L},Real{0.4L},Real{0.6L}});

  std::cout << "P % b = " << P % b << '\n';
  num_errors += (P % b) != Poly({Real{0}});

  // Pre scalar...
  Real a = 2;
  std::cout << "\na = " << a << '\n';

  std::cout << "a + Q = " << a + Q << '\n';
  num_errors += (a + Q) != Poly({Real{4},Real{1}});

  std::cout << "a - Q = " << a - Q << '\n';
  num_errors += (a - Q) != Poly({Real{0},Real{-1}});

  std::cout << "a * Q = " << a * Q << '\n';
  num_errors += (a * Q) != Poly({Real{4},Real{2}});

  std::cout << "a / Q = " << a / Q << '\n';
  num_errors += (a / Q) != Poly({0});

  std::cout << "a % Q = " << a % Q << '\n';
  num_errors += (a % Q) != Poly({2});


  Poly B;// = b;
  B = b;
  std::cout << "\nB = " << B << '\n';
  std::cout << "P % B = " << P % B << '\n';

  Q = {Real{0}, Real{-2}, Real{4}, Real{-6}, Real{8}, Real{-12}};
  std::cout << "\nQ = " << Q << '\n';

  Poly P2;
  P2 = P;
  std::cout << "\nP2 = " << P2 << '\n';
  std::cout << "P2 == P = " << (P2 == P) << '\n';

  std::cout << '\n';
  for (int i = 0; i <= 100; ++i)
    {
      Real x = i * 0.1;
      std::cout << "P(" << x << ") = " << P(x) << '\n';
    }

  CPoly
  CP({Cmplx(Real{0}, Real{-1}),
      Cmplx(Real{1}, Real{-2}),
      Cmplx(Real{2}, Real{-3}),
      Cmplx(Real{3}, Real{-4})});
  std::cout << "\nCP = " << CP << '\n';
  std::cout << "CP * CP = " << CP * CP << '\n';

  emsr::Polynomial<int> IP({0, 1, 2, 3});
  std::cout << "\nIP = " << IP << '\n';
  std::cout << "IP * IP = " << IP * IP << '\n';

  std::array<Real, 10> arr;
  P.eval(Real{1}, arr);
  std::cout << "\nP(" << Real{1} << ") =";
  for (unsigned i = 0; i < arr.size(); ++i)
    std::cout << " " << arr[i];
  std::cout << '\n';

  P.eval(Real{1}, arr.begin(), arr.end());
  std::cout << "P(" << Real{1} << ") =";
  for (auto iarr = arr.cbegin(); iarr != arr.cend(); ++iarr)
    std::cout << " " << *iarr;
  std::cout << '\n';

  std::istringstream is("(-2.0, -1.0, 0.0)");
  Poly R;
  is >> R;
  std::cout << "\nR = " << R << '\n';

  std::istringstream is2("(5.0)");
  Poly S;
  is2 >> S;
  std::cout << "\nS = " << S << '\n';

  std::istringstream is3("42.0");
  Poly T;
  is3 >> T;
  std::cout << "\nT = " << T << '\n';

  Poly u({Real{1}, Real{3}, Real{3}, Real{1}});
  Poly v({Real{1}, Real{1}});
  Poly q, r;
  emsr::divmod(u, v, q, r);
  std::cout << "\nu = " << u << '\n';
  std::cout << "v = " << v << '\n';
  std::cout << "q = " << q << '\n';
  std::cout << "r = " << r << '\n';

  Poly u1({Real{1}, Real{-3}, Real{3}, Real{-1}});
  Poly v1({Real{1}, Real{-2}, Real{1}});
  Poly q1, r1;
  emsr::divmod(u1, v1, q1, r1);
  std::cout << "\nu1 = " << u1 << '\n';
  std::cout << "v1 = " << v1 << '\n';
  std::cout << "q1 = " << q1 << '\n';
  std::cout << "r1 = " << r1 << '\n';

  Poly u2({Real{1}, Real{1}});
  Poly v2({Real{1}, Real{3}, Real{3}, Real{1}});
  Poly q2, r2;
  emsr::divmod(u2, v2, q2, r2);
  std::cout << "\nu2 = " << u2 << '\n';
  std::cout << "v2 = " << v2 << '\n';
  std::cout << "q2 = " << q2 << '\n';
  std::cout << "r2 = " << r2 << '\n';

  Poly u3({Real{1}, Real{0}, Real{0}, Real{1}});
  Poly v3({Real{1}, Real{3}, Real{3}, Real{1}});
  Poly q3, r3;
  emsr::divmod(u3, v3, q3, r3);
  std::cout << "\nu3 = " << u3 << '\n';
  std::cout << "v3 = " << v3 << '\n';
  std::cout << "q3 = " << q3 << '\n';
  std::cout << "r3 = " << r3 << '\n';

  std::cout << "P = " << P << '\n';
  std::cout << "P' = " << P.derivative() << '\n';
  std::cout << "I = " << P.integral(Real{42.0}) << '\n';

  emsr::Polynomial<Poly>
  Pp({{Real{9}, Real{8}, Real{7}}, {Real{6}, Real{5}, Real{4}}, {Real{3}, Real{2}, Real{1}}});
  std::cout << "\nPp = " << Pp << '\n';
  std::cout << "Pp(-1) = " << Pp(Real{-1.0}) << '\n';
  std::cout << "Pp(2) = " << Pp(Real{2.0}) << '\n';

  //std::array<Real, 5> aaa{{1.1, 2.2, 3.3, 4.4, 5.5}};
  //std::cout << "aaa = " << emsr::Polynomial_eval(aaa, 3.13149) << '\n';

  std::cout << "\nP = " << P << '\n';
  std::cout << "P(1) = " << P(Real{1.0}) << '\n';
  // This can fail compile if complex literals return something other than std::complex.
  //if (std::is_same_v<decltype(1.0i), Cmplx>)
  //  std::cout << "P(i) = " << P(1.0i) << '\n';
  std::cout << "P(i) = " << P(Cmplx{0, 1}) << '\n';

  Poly e([](unsigned int k)
         -> Real
         { return 1.0 / factorial<Real>(k); }, 20);
  std::cout << "\ne = " << e << '\n';
  std::cout << "e(1) = " << e(1) << '\n';
  std::cout << "e(1) - exp(1) = " << e(1) - std::exp(Real(1)) << '\n';

  return num_errors;
}

