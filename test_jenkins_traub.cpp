/*
$HOME/bin/bin/g++ -std=gnu++20 -g -Wall -Wextra -Wno-psabi -I. -I../include -o test_jenkins_traub test_jenkins_traub.cpp -lquadmath
*/

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "ext/polynomial.h"
#include "ext/solver_jenkins_traub.h"


template<typename _Real>
  void
  test_jenkins_traub()
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);

    int order = 0;
    int MAX_TERMS = 1000;
    while ((order < 2) || (order > MAX_TERMS - 1))
      {
	std::cout << "Polynomial order (2 - " << MAX_TERMS - 1 << "): ";
	std::cin >> order;
      }
    std::vector<_Real> a(order + 1);

    std::cout << "Enter coefficients, high order to low order.\n";
    for (int i = 0; i <= order; ++i)
      {
	std::cout << "a[" << (order - i) << "] = ";
	std::cin >> a[i];
      }
    std::cout << "\nP : ( ";
    for (int i = 0; i < order; ++i)
      std::cout << a[order - i] << ", ";
    std::cout << a[0] << " )\n";

    __gnu_cxx::_JenkinsTraubSolver jenkins_traub(a);
    const auto zeros = jenkins_traub.solve();
    std::cout << "\nThe zeros are:\n";
    for (const auto& z : zeros)
      std::cout << z << '\n';

    std::cout << "\nSolution tests:\n";
    // Remember to reverse the polynomial coefficients!
    __gnu_cxx::_Polynomial<_Real> poly(a.rbegin(), a.rend());
    for (const auto& z : zeros)
      {
	const auto idx = z.index();
	std::cout << "P(" << z << ") = ";
	if (idx == 1)
	  std::cout << poly(std::get<1>(z));
	else if (idx == 2)
	  std::cout << poly(std::get<2>(z));
	else
	  std::cout << "error";
	std::cout << '\n';
      }
  }

int
main()
{
  std::cout << "\ndouble\n======\n" << std::flush;
  test_jenkins_traub<double>();

  std::cout << "\nlong double\n===========\n" << std::flush;
  test_jenkins_traub<long double>();

  std::cout << "\nfloat\n=====\n" << std::flush;
  test_jenkins_traub<float>();
}

