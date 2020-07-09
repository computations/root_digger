#ifndef RD_TEST_UTIL_HPP
#define RD_TEST_UTIL_HPP
#include <string>

constexpr char base_58_chars[] =
    "123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz";

inline size_t compute_digit_with_base(size_t i, size_t n, size_t base);

std::string base_58_encode(uint32_t n);

#endif
