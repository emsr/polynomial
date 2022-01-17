
#include <iostream>
#include <iomanip>

#include <ext/polynomial.h>
#include <ext/solver_jenkins_traub.h>


template<typename Real>
  int
  test_deriv_integ()
  {
    int num_errors = 0;

    emsr::Polynomial<Real> poly({5, -4, 3, -2, 1});
    const auto poly_deriv = poly.derivative();
    const auto poly_integ = poly.integral(-6);

    for (int i = -2; i <= 2; ++i)
    {
        const auto c = Real(i);
        const auto dp1 = poly_deriv(c);
        const auto dp2 = poly.derivative(c);
        num_errors += dp1 != dp2;

        const auto a = Real(i);
        const auto b = Real(2*i);
        const auto ip1 = poly_integ(b) - poly_integ(a);
        const auto ip2 = poly.integral(a, b);
        num_errors += ip1 != ip2;
    }

    return num_errors;
  }

int
main()
{
  int num_errors = 0;

  std::cout << "\ndouble\n======\n" << std::flush;
  num_errors += test_deriv_integ<double>();

  std::cout << "\nlong double\n===========\n" << std::flush;
  num_errors += test_deriv_integ<long double>();

  std::cout << "\nfloat\n=====\n" << std::flush;
  num_errors += test_deriv_integ<float>();

  return num_errors;
}

