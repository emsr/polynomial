/*
$HOME/bin/bin/g++ -std=gnu++20 -g -Wall -Wextra -Wno-psabi -I. -I../include -o test_solution test_solution.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_solution > test_solution.txt
*/

#include <ext/solution.h>

int
main()
{
  __gnu_cxx::solution_t<double> s1(3.0);
  auto s2 = __gnu_cxx::solution_t<double>(std::complex<double>{1.0, 3.0});
  auto x1 = 4.0 * s1;
}
