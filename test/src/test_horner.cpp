
#include "ext/horner.h"

#include <iostream>
#include <iomanip>

int
main()
{
  std::cout << __gnu_cxx::horner_big_end(3.0, 1.0, 2.0, 1.0) << '\n';
  std::cout << __gnu_cxx::horner_big_end(3.0, 1.0, 2.0, 3.0) << '\n';
  std::cout << __gnu_cxx::horner(3.0, 3.0, 2.0, 1.0) << '\n';
}
