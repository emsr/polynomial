
#include <limits>
#include <complex>
#include <iostream>

#include <emsr/polynomial.h>
#include <emsr/solver_laguerre.h>

template<typename Tp>
  void
  test_laguerre_step()
  {
    using cmplx = std::complex<Tp>;

    emsr::Polynomial<cmplx>
    CP({cmplx(0.0, -1.0), cmplx(1.0, -2.0),
	cmplx(2.0, -3.0), cmplx(3.0, -4.0)});

    std::cout << "\ntest_laguerre - stepper";
    std::cout << "\n-----------------------\n";
    std::cout << "CP    = " << CP << '\n';
    auto lagsolver =  emsr::LaguerreSolver(CP);
    //auto roots = lagsolver.solve();
    //for (x : roots)
    for (unsigned int i = 0; i < CP.degree(); ++i)
      {
	std::cout << '\n';
	auto x = lagsolver.step();
	std::cout << "x     = " << x << '\n';
	std::cout << "iters = " << lagsolver.num_iters() << '\n';
	std::cout << "CP(x) = " << CP(x) << '\n';
	std::cout << "CP    = " << lagsolver.polynomial() << '\n';
      }

    std::cout << "\ntest_laguerre - full solver";
    std::cout << "\n---------------------------\n";
    auto lagsolver2 =  emsr::LaguerreSolver(CP);
    auto zeros = lagsolver2.solve();
    for(auto& z : zeros)
      std::cout << "z     = " << z << '\n';
  }

int
main()
{
  test_laguerre_step<double>();
}
