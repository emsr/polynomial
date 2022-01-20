
#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/polynomial.h>
#include <emsr/solver_jenkins_traub.h>


template<typename Real>
  int
  test_jenkins_traub()
  {
    std::cout.precision(std::numeric_limits<Real>::digits10);

    int order = 0;
    int MAX_TERMS = 1000;
    while ((order < 2) || (order > MAX_TERMS - 1))
      {
	std::cout << "Polynomial order (2 - " << MAX_TERMS - 1 << "; 0 to quit): ";
	std::cin >> order;
        if (order <= 0)
          break;
      }
    std::vector<Real> a(order + 1);

    std::cout << "Enter coefficients, high order to low order.\n";
    for (int i = 0; i <= order; ++i)
      {
	std::cout << "a[" << (order - i) << "] = ";
	std::cin >> a[i];
      }
    std::cout << "\nP : ( ";
    for (int i = 0; i <= order; ++i)
      std::cout << a[order - i] << ", ";
    std::cout << a[0] << " )\n";

    emsr::JenkinsTraubSolver jenkins_traub(a);
    const auto zeros = jenkins_traub.solve();
    std::cout << "\nThe zeros are:\n";
    for (const auto& z : zeros)
      std::cout << z << '\n';

    int num_errors = 0;
    const auto tol = std::sqrt(std::numeric_limits<Real>::epsilon());

    std::cout << "\nSolution tests:\n";
    // Remember to reverse the polynomial coefficients!
    emsr::Polynomial<Real> poly(a.rbegin(), a.rend());
    for (const auto& z : zeros)
      {
	const auto idx = z.index();
	std::cout << "P(" << z << ") = ";
	if (idx == 1)
	  {
	    const auto P = poly(std::get<1>(z));
	    std::cout << P << '\n';
	    if (const auto aP = std::abs(P); aP > tol)
	      {
		++num_errors;
		std::cout << "Fail: |P| = " << aP << '\n';
	      }
	  }
	else if (idx == 2)
	  {
	    const auto P = poly(std::get<2>(z));
	    std::cout << P << '\n';
	    if (const auto aP = std::abs(P); aP > tol)
	      {
		++num_errors;
		std::cout << "Fail: |P| = " << aP << '\n';
	      }
	  }
	else
	  {
	    std::cout << "error\n";
	  }
      }

    return num_errors;
  }

int
main()
{
  int num_errors = 0;

  std::cout << "\ndouble\n======\n" << std::flush;
  num_errors += test_jenkins_traub<double>();

  std::cout << "\nlong double\n===========\n" << std::flush;
  num_errors += test_jenkins_traub<long double>();

  std::cout << "\nfloat\n=====\n" << std::flush;
  num_errors += test_jenkins_traub<float>();

  return num_errors;
}

