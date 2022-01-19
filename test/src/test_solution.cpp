
#include <emsr/solution.h>

int
main()
{
  emsr::solution_t<double> s1(3.0);
  auto s2 = emsr::solution_t<double>(std::complex<double>{1.0, 3.0});
  [[maybe_unused]] auto x1 = 4.0 * s1;
  [[maybe_unused]] auto x2 = 4.0 * s2;

  [[maybe_unused]] auto b1 = (s1 == 3.0);
}
