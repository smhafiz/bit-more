#include <iostream>
#include "constexpr_for.h"

namespace {

void double_loop() {
  constexpr size_t array_size = 3;
  std::array<std::array<size_t, array_size>, array_size> zero{
      {{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}};
    auto values = zero;
    for_constexpr<for_bounds<0, 3>, for_bounds<0, 4>>(
      [&](auto i, auto j) 
      { //values[i][j]++;
        std::cout<< i << ", " << j << ", " << values[i][j] << std::endl;
      }
    );

}
} // namespace

int main() {
  double_loop();
}