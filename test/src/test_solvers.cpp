
#include <experimental/array>
#include <iostream>
#include <iomanip>
#include <limits>

#include <emsr/solver_low_degree.h>

template<typename Real>
  struct status
  {
    bool ok;
    Real err;
  };

template<typename Real, unsigned long int Dim>
  status<Real>
  try_solution(std::array<Real, Dim> coef,
	       const emsr::Solution<Real>& z)
  {
    const auto s_eps = Real{100} * std::numeric_limits<Real>::epsilon();
    if (z.index() == 0)
      return {true, Real{0}};
    else if (z.index() == 1)
      {
	const Real x = std::get<1>(z);
	Real y = coef[Dim - 1];
	for (unsigned int i = Dim - 1; i > 0; --i)
	  y = coef[i - 1] + x * y;
	const Real eps = std::abs(y);
	return {std::abs(y) < s_eps * std::abs(x), eps};
      }
    else if (z.index() == 2)
      {
	const std::complex<Real> x = std::get<2>(z);
	std::complex<Real> y = coef[Dim - 1];
	for (unsigned int i = Dim - 1; i > 0; --i)
	  y = coef[i - 1] + x * y;
	const Real eps = std::abs(y);
	return {std::abs(y) < s_eps * std::abs(x), eps};
      }
    return {false, Real{0}};
  }

int
main()
{
  using std::experimental::make_array;

  std::cout << std::boolalpha;
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto w = 8 + std::cout.precision();

  // 0 = (z + 2)(z - 1) = z^2 + z - 2
  std::cout << '\n';
  const auto quad_coef = make_array(-2.0, 1.0, 1.0);
  auto quad = emsr::quadratic<double>(quad_coef);
  for (int i = 0; i < 2; ++i)
    {
      auto stat = try_solution(quad_coef, quad[i]);
      std::cout << std::setw(w) << quad[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  // 0 = [z - (2+i)][z - (2-i)] = z^2 -4z + 5
  std::cout << '\n';
  const auto cquad_coef = make_array(5.0, -4.0, 1.0);
  auto cquad = emsr::quadratic<double>(cquad_coef);
  for (int i = 0; i < 2; ++i)
    {
      auto stat = try_solution(cquad_coef, cquad[i]);
      std::cout << std::setw(w) << cquad[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  // 0 = z - 2 = 0z^2 + z - 2
  std::cout << '\n';
  const auto linquad_coef = make_array(-2.0, 1.0, 0.0);
  auto linquad = emsr::quadratic<double>(linquad_coef);
  for (int i = 0; i < 2; ++i)
    {
      auto stat = try_solution(linquad_coef, linquad[i]);
      std::cout << std::setw(w) << linquad[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  // 0 = (2z + 6)(z - 4)(z - 1) = 2z^3 - 4z^2 - 22z + 24
  std::cout << '\n';
  const auto cub_coef = make_array(24.0, -22.0, -4.0, 2.0);
  auto cub = emsr::cubic<double>(cub_coef);
  for (int i = 0; i < 3; ++i)
    {
      auto stat = try_solution(cub_coef, cub[i]);
      std::cout << std::setw(w) << cub[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  std::cout << '\n';
  const auto quadcub_coef = make_array(24.0, -22.0, -4.0, 0.0);
  auto quadcub = emsr::cubic<double>(quadcub_coef);
  for (int i = 0; i < 3; ++i)
    {
      auto stat = try_solution(quadcub_coef, quadcub[i]);
      std::cout << std::setw(w) << quadcub[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  std::cout << '\n';
  const auto lincub_coef = make_array(24.0, -22.0, 0.0, 0.0);
  auto lincub = emsr::cubic<double>(lincub_coef);
  for (int i = 0; i < 3; ++i)
    {
      auto stat = try_solution(lincub_coef, lincub[i]);
      std::cout << std::setw(w) << lincub[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  std::cout << '\n';
  const auto cub2_coef = make_array(-4.0, -15.0, 0.0, 1.0);
  auto cub2 = emsr::cubic<double>(cub2_coef);
  for (int i = 0; i < 3; ++i)
    {
      auto stat = try_solution(cub2_coef, cub2[i]);
      std::cout << std::setw(w) << cub2[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  // 0 = (z - 5)[z - (2+i)][z - (2-i)] = z^3 - 9z^2 + 25z - 25
  std::cout << '\n';
  const auto ccube_coef = make_array(-25.0, 25.0, -9.0, 1.0);
  auto ccube = emsr::cubic<double>(ccube_coef);
  for (int i = 0; i < 3; ++i)
    {
      auto stat = try_solution(ccube_coef, ccube[i]);
      std::cout << std::setw(w) << ccube[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  std::cout << '\n';
  const auto quart_coef = make_array(-30.0, 31.0, 5.0, -7.0, 1.0);
  auto quart = emsr::quartic<double>(quart_coef);
  for (int i = 0; i < 4; ++i)
    {
      auto stat = try_solution(quart_coef, quart[i]);
      std::cout << std::setw(w) << quart[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  std::cout << '\n';
  const auto cubquart_coef = make_array(-30.0, 31.0, 5.0, -7.0, 0.0);
  auto cubquart = emsr::quartic<double>(cubquart_coef);
  for (int i = 0; i < 4; ++i)
    {
      auto stat = try_solution(cubquart_coef, cubquart[i]);
      std::cout << std::setw(w) << cubquart[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  std::cout << '\n';
  const auto quadquart_coef = make_array(-30.0, 31.0, 5.0, 0.0, 0.0);
  auto quadquart = emsr::quartic<double>(quadquart_coef);
  for (int i = 0; i < 4; ++i)
    {
      auto stat = try_solution(quadquart_coef, quadquart[i]);
      std::cout << std::setw(w) << quadquart[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  std::cout << '\n';
  const auto linquart_coef = make_array(-30.0, 31.0, 0.0, 0.0, 0.0);
  auto linquart = emsr::quartic<double>(linquart_coef);
  for (int i = 0; i < 4; ++i)
    {
      auto stat = try_solution(linquart_coef, linquart[i]);
      std::cout << std::setw(w) << linquart[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  // 0 = (z - 3)(z - 5)[z - (2+i)][z - (2-i)]
  //   = z^4 - 12z^3 + 52z^2 - 100z + 75
  std::cout << '\n';
  const auto cquart_coef = make_array(75.0, -100.0, 52.0, -12.0, 1.0);
  auto cquart = emsr::quartic<double>(cquart_coef);
  for (int i = 0; i < 4; ++i)
    {
      auto stat = try_solution(cquart_coef, cquart[i]);
      std::cout << std::setw(w) << cquart[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  // 0 = [z - (-1+i)][z - (-1-i)][z - (4+2i)][z - (4-2i)]
  //   = z^4 - 6z^3 + 6z^2 + 24z + 40
  std::cout << '\n';
  const auto biquad_coef = make_array(40.0, 24.0, 6.0, -6.0, 1.0);
  auto biquad = emsr::quartic<double>(biquad_coef);
  for (int i = 0; i < 4; ++i)
    {
      auto stat = try_solution(biquad_coef, biquad[i]);
      std::cout << std::setw(w) << biquad[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  // 0 = (z - 1)(z + 1)(z - 4)(z + 4)
  //   = (z^2 - 1)(z^2 - 16)
  //   = z^4 - 17z^2 + 16
  std::cout << '\n';
  const auto biquad2_coef = make_array(16.0, 0.0, -17.0, 0.0, 1.0);
  auto biquad2 = emsr::quartic<double>(biquad2_coef);
  for (int i = 0; i < 4; ++i)
    {
      auto stat = try_solution(biquad2_coef, biquad2[i]);
      std::cout << std::setw(w) << biquad2[i] << "  good: " << std::setw(5) << stat.ok << "  error: " << stat.err << '\n';
    }

  // New API.
  emsr::quadratic<double>(-2.0, 1.0, 1.0);
  emsr::cubic<double>(24.0, -22.0, -4.0, 2.0);
  emsr::quartic<double>(-30.0, 31.0, 5.0, -7.0, 1.0);
}
