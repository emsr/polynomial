
#include <ext/horner.h>

#include <iostream>
#include <iomanip>

int
main()
{
  int num_errors = 0;

  const auto p1_3 = __gnu_cxx::horner_big_end(3.0, 1.0, 2.0, 1.0);
  num_errors += p1_3 != 16.0;
  std::cout << "x^2 + 2x + 1; x = 3: " << p1_3 << '\n';

  const auto p2_3 = __gnu_cxx::horner_big_end(3.0, 1.0, 2.0, 3.0);
  num_errors += p2_3 != 18.0;
  std::cout << "x^2 + 2x + 3; x = 3: " << p2_3 << '\n';

  const auto p3_3 = __gnu_cxx::horner(3.0, 3.0, 2.0, 1.0);
  num_errors += p3_3 != 18.0;
  std::cout << "3 + 2x + x^2; x = 3: " << p3_3 << '\n';

  return num_errors;
}
