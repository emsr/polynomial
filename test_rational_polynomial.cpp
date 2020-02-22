
#include <iostream>
#include <complex>
#include <sstream>

#include "ext/rational_polynomial.h"

int
main()
{
  __gnu_cxx::_Polynomial<double> P({0.0, 1.0, 2.0, 3.0});
  __gnu_cxx::_Polynomial<double> Q({2.0, 1.0});
  __gnu_cxx::_RationalPolynomial<double> R(P, Q);

  std::cout << "R = " << R << '\n';

  std::cout << "P = " << R.numer() << '\n';
  std::cout << "+P = " << +R.numer() << '\n';
  std::cout << "-P = " << -R.numer() << '\n';
  std::cout << "P = " << R.numer() << '\n';
  std::cout << "degree(P) = " << R.numer().degree() << '\n';

  std::cout << "Q = " << R.denom() << '\n';
  std::cout << "+Q = " << +R.denom() << '\n';
  std::cout << "-Q = " << -R.denom() << '\n';
  std::cout << "Q = " << R.denom() << '\n';
  std::cout << "degree(Q) = " << R.denom().degree() << '\n';

}

