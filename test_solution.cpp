
#include <ext/solution.h>

int
main()
{
  __gnu_cxx::solution_t<double> s1(3.0);
  auto s2 = __gnu_cxx::solution_t<double>(std::complex<double>{1.0, 3.0});
  auto x1 = 4.0 * s1;
}
