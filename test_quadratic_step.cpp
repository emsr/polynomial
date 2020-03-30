
#include <limits>
#include <complex>
#include <ext/polynomial.h>
#include <ext/solver_quadratic.h>

#include <iostream>

template<typename _Tp>
  void
  test_quadratic_step()
  {
    __gnu_cxx::_Polynomial<std::complex<_Tp>>
    CP({std::complex<_Tp>(0.0, -1.0),
	std::complex<_Tp>(1.0, -2.0),
	std::complex<_Tp>(2.0, -3.0),
	std::complex<_Tp>(3.0, -4.0)});

    std::cout << "\ntest_quadratic\n";
    std::cout << "CP = " << CP << '\n';

    __gnu_cxx::_QuadraticSolver<_Tp> qsolve(CP);
    for (unsigned i = 0; i < CP.degree(); i += 2)
      {
	std::cout << '\n';
	std::complex<_Tp> b, c;
	auto q = qsolve.step();
	std::cout << "b   = " << q[1] << '\n';
	std::cout << "c   = " << q[0] << '\n';
	std::cout << "CP = " << qsolve.polynomial() << '\n';
      }
  }

int
main()
{
  test_quadratic_step<double>();
}
