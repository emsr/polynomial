/*
$HOME/bin/bin/g++ -std=c++2a -Iinclude -g -Wall -Wextra -o test_laguerre_step test_laguerre_step.cpp
*/

#include <limits>
#include <complex>
#include <iostream>

#include <ext/polynomial.h>
#include <ext/solver_laguerre.h>

template<typename _Tp>
  void
  test_laguerre_step()
  {
    using cmplx = std::complex<_Tp>;

    __gnu_cxx::_Polynomial<cmplx>
    CP({cmplx(0.0, -1.0), cmplx(1.0, -2.0),
	cmplx(2.0, -3.0), cmplx(3.0, -4.0)});

    std::cout << "\ntest_laguerre - stepper";
    std::cout << "\n-----------------------\n";
    std::cout << "CP    = " << CP << '\n';
    auto lagsolver =  __gnu_cxx::_LaguerreSolver(CP);
    //auto roots = lagsolver.solve();
    //for (x : roots)
    for (unsigned int i = 0; i < CP.degree(); ++i)
      {
	std::cout << '\n';
	auto x = lagsolver.step();
	std::cout << "x     = " << x << '\n';
	std::cout << "its   = " << lagsolver.num_iters() << '\n';
	std::cout << "CP(x) = " << CP(x) << '\n';
	std::cout << "CP    = " << lagsolver.polynomial() << '\n';
      }

    std::cout << "\ntest_laguerre - full solver";
    std::cout << "\n---------------------------\n";
    auto lagsolver2 =  __gnu_cxx::_LaguerreSolver(CP);
    auto zeros = lagsolver2.solve();
    for(auto& z : zeros)
      std::cout << "z     = " << z << '\n';
  }

int
main()
{
  test_laguerre_step<double>();
}
