
#include <ext/polynomial.h>

#include <iostream>
#include <iomanip>

int
main()
{
  emsr::Polynomial<double> poly({-1, 0, 0, 0, 0, 1});
  std::cout << poly << '\n';
  poly.shift(2);
  std::cout << poly << '\n';
  emsr::Polynomial<double> new_poly({31, 80, 80, 40, 10, 1});
  std::cout << new_poly << '\n';
  return poly != new_poly;
}
