#include <cte>
#include <iostream>
using namespace cte_random;
int main() {

  int s = 1234;

  std::printf("%12.6e\n", ran3(&s));
  return 0;
}