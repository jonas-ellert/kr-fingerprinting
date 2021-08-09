#pragma once

#include "kr-fingerprinting64.hpp"

namespace kr_fingerprinting {

namespace u128 {

constexpr uint64_t countr_one(uint128_t v) {
  auto r = std::countr_one((uint64_t)v);
  return (r < 64) ? r : (std::countr_one((uint64_t)(v >> 64)) + 64);
}

constexpr uint64_t popcount(uint128_t v) {
  return std::popcount((uint64_t)v) + std::popcount((uint64_t)(v >> 64));
}

template <uint128_t modulus>
inline constexpr uint128_t mod(uint128_t const value) {
  static_assert(popcount(modulus) < 128);
  static_assert(popcount(modulus) == countr_one(modulus));
  constexpr uint64_t shift = popcount(modulus);
  uint128_t const i = (value & modulus) + (value >> shift);
  return (i >= modulus) ? (i - modulus) : i;
}

// this assumes a,b,c < modulus
// a * b + c
template <uint128_t modulus>
inline constexpr static uint128_t mult_add(uint128_t const a, uint128_t const b,
                                           uint128_t const c) {
  if (a >= modulus) __builtin_unreachable();
  if (b >= modulus) __builtin_unreachable();
  if (c >= modulus) __builtin_unreachable();

  constexpr uint64_t shift = popcount(modulus);

  uint128_t const al = (uint64_t)a;
  uint128_t const ah = a >> 64;
  uint128_t const bl = (uint64_t)b;
  uint128_t const bh = b >> 64;

  uint128_t const h = ah * bh;
  uint128_t const m1 = ah * bl;
  uint128_t const m2 = bh * al;
  uint128_t const l = al * bl;

  if constexpr (shift < 127) {
    if (h >= modulus) __builtin_unreachable();
    if (l >= modulus) __builtin_unreachable();
    if (m1 >= modulus) __builtin_unreachable();
    if (m2 >= modulus) __builtin_unreachable();
    // this only works because we have sufficiently many overflow bits
    uint128_t const m = m1 + m2;
    uint128_t const sum = c + (l & modulus) + (l >> shift) +
                          ((m << 64) & modulus) + (m >> (shift - 64)) +
                          ((h << (128 - shift)) & modulus) +
                          (h >> (2 * shift - 128));
    return mod<modulus>(sum);
  } else {
    // for p127 the single overflow bit is not sufficient
    uint128_t const carry =
        ((l >> 64) + (c >> 64) + (uint64_t)m1 + (uint64_t)m2) >> 64;

    uint128_t const h128 = h + (m1 >> 64) + (m2 >> 64) + carry;
    uint128_t const l128 = l + c + (m1 << 64) + (m2 << 64);

    uint128_t sum =
        ((h128 << (128 - shift)) | (l128 >> shift)) + (l128 & modulus);

    return mod<modulus>(sum);
  }
}

template <uint128_t modulus>
inline constexpr static uint128_t mult(uint128_t const a, uint128_t const b) {
  return mult_add<modulus>(a, b, 0);
}

// fast squaring
template <uint128_t modulus>
inline constexpr uint128_t power(uint128_t base, uint128_t exponent) {
  uint128_t result = 1;
  uint128_t b = base;
  while (exponent > 0) {
    if (exponent & 1ULL) result = mult<modulus>(b, result);
    b = mult<modulus>(b, b);
    exponent >>= 1;
  }
  return result;
}

inline static uint128_t random(uint128_t min, uint128_t max) {
  constexpr uint64_t max64 = std::numeric_limits<uint64_t>::max();

  max -= min;
  if (max <= max64) return min + u64::random(0, max);

  uint64_t const shift = 64 - std::bit_width((uint64_t)(max >> 64));
  uint64_t const a = u64::random(0, max64);
  uint64_t const b = u64::random(0, max64);
  uint128_t const r = (((((uint128_t)a) << 64) | b) << shift) >> shift;

  return min + ((r > max) ? random(0, max) : r);
}

template <uint64_t shift>
struct sliding_windowX {
 private:
  // mersenne prime 2^61 - 1
  constexpr static uint64_t s = shift;
  constexpr static uint128_t p = (((uint128_t)1) << s) - 1;

  uint64_t const window_size_;
  uint128_t const base_;
  double const collision_rate_ = ((double)window_size_ - 1) / p;

  table2d<uint128_t> const table_ = table2d<uint128_t>([&](uint128_t d[][256]) {
    uint128_t const max_exponent = u128::power<p>(base_, window_size_);
    for (uint64_t i = 0; i < 256; ++i) {
      d[i][0] = u128::mod<p>(p - u128::mult<p>(i, max_exponent));
      for (uint64_t j = 1; j < 256; ++j) {
        d[i][j] = u128::mod<p>(d[i][j - 1] + 1);
      }
    }
  });

 public:
  sliding_windowX(uint64_t const window_size, uint128_t const base)
      : window_size_(window_size), base_(u128::mod<p>(base)){};

  sliding_windowX(uint64_t const window_size)
      : sliding_windowX(window_size, u128::random(1, p - 1)){};

  template <ByteType T>
  inline uint128_t roll_right(uint128_t const fp, T const pop_left,
                              T const push_right) const {
    auto lookup = table_(pop_left, push_right);
    if (base_ >= p || fp >= p || lookup >= p)
      __builtin_unreachable();
    else
      return u128::mult_add<p>(base_, fp, lookup);
  }

  template <ByteType T>
  inline uint128_t roll_right(uint128_t const fp, T const push_right) const {
    // if (base_ >= p || fp >= p)
    //  __builtin_unreachable();
    // else
    return u128::mult_add<p>(base_, fp, push_right);
  }

  inline uint128_t base() const { return base_; }
  inline uint64_t window_size() const { return window_size_; }
  inline uint64_t bits() const { return s; }
  inline double collision_rate() const { return collision_rate_; }
};

}  // namespace u128

}  // namespace kr_fingerprinting
