
#include <limits>
#include <complex>
#include <iostream>

#include <emsr/polynomial.h>
#include <emsr/solver_quadratic.h>

template<typename Tp>
  void
  test_quadratic_step()
  {
    emsr::Polynomial<std::complex<Tp>>
    CP({std::complex<Tp>(0.0, -1.0),
	std::complex<Tp>(1.0, -2.0),
	std::complex<Tp>(2.0, -3.0),
	std::complex<Tp>(3.0, -4.0)});

    std::cout << "\ntest_quadratic\n";
    std::cout << "CP = " << CP << '\n';

    emsr::QuadraticSolver<Tp> qsolve(CP);
    for (unsigned i = 0; i < CP.degree(); i += 2)
      {
	std::cout << '\n';
	std::complex<Tp> b, c;
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
