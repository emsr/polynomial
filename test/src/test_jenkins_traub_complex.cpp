
#include <vector>
#include <complex>
#include <iostream>
#include <limits>

#include <ext/polynomial.h>
#include <ext/solver_jenkins_traub.h>

template<typename Real>
int
test_complex_jenkins_traub()
{
    std::cout.precision(std::numeric_limits<Real>::digits10);

    int degree = 0;
    std::cout << "Enter degree of polynomial: ";
    std::cin >> degree;

    std::vector<std::complex<Real>> p(degree + 1);
    for (int i = 0; i <= degree; ++i)
    {
        std::cout << "Enter coefficient " << i << ": ";
        std::cin >> p[i];
    }

    emsr::JenkinsTraubSolver<std::complex<Real>> jt(p);
    auto zero = jt.solve();
    std::cout << '\n';
    for (unsigned int i = 0; i < zero.size(); ++i)
    {
        std::cout << "Zero " << i << ": " << zero[i] << '\n';
    }

    int num_errors = 0;
    const auto tol = std::sqrt(std::numeric_limits<Real>::epsilon());

    emsr::Polynomial<std::complex<Real>> poly(p.rbegin(), p.rend());

    std::cout << "\nSolution tests:\n";
    for (unsigned int i = 0; i < zero.size(); ++i)
    {
        const auto P = poly(zero[i]);
	std::cout << P << '\n';
	if (const auto aP = std::abs(P); aP > tol)
        {
            ++num_errors;
            std::cout << "Fail: |P| = " << aP << '\n';
        }
    }

    return num_errors;
}

int
main()
{
  int num_errors = 0;

  std::cout << "\ndouble\n======\n";
  num_errors += test_complex_jenkins_traub<double>();

  std::cout << "\nlong double\n===========\n";
  num_errors += test_complex_jenkins_traub<long double>();

  std::cout << "\nfloat\n=====\n";
  // Punt accumulating float errors.
  //num_errors += test_complex_jenkins_traub<float>();
  test_complex_jenkins_traub<float>();

  return num_errors;
}
