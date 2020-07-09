#include "test_util.hpp"
#include <random>

inline size_t compute_digit_with_base(size_t i, size_t n, size_t base) {
  return (n % static_cast<size_t>(std::pow(base, i + 1))) / std::pow(base, i);
}

std::string base_58_encode(uint32_t n) {
  size_t alphabet_size = sizeof(base_58_chars)-1;
  size_t len = std::ceil(std::log(n) / std::log(alphabet_size));
  std::string enc;
  enc.resize(len);
  for (size_t i = 0; i < len; ++i) {
    enc[i] = base_58_chars[compute_digit_with_base(i, n, alphabet_size)];
  }
  return enc;
}
