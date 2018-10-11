/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -I.. -o test_jenkins_traub test_jenkins_traub.cpp -lquadmath
*/

#include <vector>
#include <limits>
#include <cmath>

#include "ext/polynomial.h"
#include "ext/solver_jenkins_traub.h"


template<typename _Real>
  void
  test_jenkins_traub()
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);

    int order = 0;
    int MAX_TERMS = 100;
    while ((order < 2) || (order > MAX_TERMS - 1))
      {
	std::cout << "Polynomial order (2 - " << MAX_TERMS - 1 << "): ";
	std::cin >> order;
      }
    std::vector<_Real> a(order + 1);

    std::cout << "Enter coefficients, high order to low order.\n";
    for (int i = 0; i <= order; ++i)
      {
	std::cout << "a[" << i << "] = ";
	std::cin >> a[i];
      }

    __gnu_cxx::_JenkinsTraubSolver jenkins_traub(a);
    const auto zeros = jenkins_traub.solve();
    std::cout << "\nThe zeros are:\n";
    for (const auto& z : zeros)
      std::cout << z << '\n';
/*
    const auto eq = jenkins_traub.equations();
    std::cout << "\nThe quadratic factors are:\n";
    for (int p = 0; p < eq.size() / 2; ++p)
      std::cout << "t^2 + " << eq[2 * p + 1] << " t + " << eq[2 * p] << '\n';
    if ((eq.size() % 2) == 1)
      std::cout << "The linear term is: \nt - " << eq.back() << '\n';
*/
    std::cout << "\nSolution tests:\n";
    __gnu_cxx::_Polynomial<_Real> poly(a.begin(), a.end());
    for (const auto& z : zeros)
      {
	const auto idx = z.index();
	std::cout << "P(" << z << ") = ";
	if (idx == 1)
	  std::cout << poly(std::get<1>(z));
	else if (idx == 2)
	  std::cout << poly(std::get<2>(z));
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

