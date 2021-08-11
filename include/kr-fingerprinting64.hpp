#pragma once

#include <bit>
#include <cstring>
#include <limits>
#include <random>
#include <sstream>
#include <type_traits>
#include "tuple.hpp"

#define KRINLNFN __attribute__((always_inline)) inline

namespace kr_fingerprinting {

__extension__ using uint128_t = unsigned __int128;

template <typename T>
concept ByteType = std::unsigned_integral<T> && (sizeof(T) == 1);

template <typename T>
inline T *const_auto_cast(T const *t) {
  return const_cast<T *>(t);
}

namespace u64 {

constexpr uint64_t p61 = (1ULL << 61) - 1;

KRINLNFN constexpr uint64_t mod(uint128_t const value) {
  uint64_t const i = (value & p61) + (value >> 61);
  return (i & p61) + (i >> 61);
}

// fast squaring
constexpr uint64_t power(uint64_t base, uint64_t exponent) {
  uint64_t result = 1;
  uint128_t b = base;
  while (exponent > 0) {
    if (exponent & 1ULL) result = mod(b * result);
    b = mod(b * b);
    exponent >>= 1;
  }
  return result;
}

inline static uint64_t random(uint64_t min, uint64_t max) {
  static std::mt19937_64 g = std::mt19937_64(std::random_device()());
  return (std::uniform_int_distribution<uint64_t>(min, max))(g);
}

struct sliding_window61 {
 private:
  uint64_t const window_size_;
  uint64_t const base_;

  uint64_t const table_[256][256] = {};

  double const collision_rate_ = ((double)window_size_ - 1) / (p61 - 2);

 public:
  using fingerprint_type = uint64_t;

  sliding_window61(uint64_t const window_size, uint64_t const base)
      : window_size_(window_size), base_(u64::mod(base)) {
    auto d = const_auto_cast(table_);
    uint64_t const max_exponent = u64::power(base_, window_size_);
    for (uint64_t i = 0; i < 256; ++i) {
      d[i][0] = u64::mod(p61 - u64::mod(i * (uint128_t)max_exponent));
      for (uint64_t j = 1; j < 256; ++j) {
        d[i][j] = u64::mod(d[i][j - 1] + 1);
      }
    }
  };

  sliding_window61(uint64_t const window_size)
      : sliding_window61(window_size, u64::random(1, p61 - 1)){};

  template <ByteType T>
  KRINLNFN uint64_t roll_right(uint64_t const fp, T const pop_left,
                               T const push_right) const {
    auto lookup = table_[pop_left][push_right];
    if (base_ >= p61 || fp >= p61 || lookup >= p61)
      __builtin_unreachable();
    else
      return u64::mod(((uint128_t)base_) * fp + lookup);
  }

  template <ByteType T>
  KRINLNFN uint64_t roll_right(uint64_t const fp, T const push_right) const {
    if (base_ >= p61 || fp >= p61)
      __builtin_unreachable();
    else
      return u64::mod(((uint128_t)base_) * fp + push_right);
  }

  inline uint64_t base() const { return base_; }
  inline uint64_t window_size() const { return window_size_; }
  inline uint64_t bits() const { return 61; }
  inline double collision_rate() const { return collision_rate_; }
};

template <uint64_t x>
struct sliding_window_multi61 {
 private:
  using tuple = kr_tuple::tuple<x>;

  uint64_t const window_size_;
  tuple const base_;

  tuple const table_[256][256] = {};

  double const collision_rate_ = std::pow(((double)window_size_ - 1) / p61, x);

 public:
  using fingerprint_type = tuple;

  sliding_window_multi61(uint64_t const window_size, tuple base)
      : window_size_(window_size), base_(base.apply(u64::mod)) {
    static_assert(x > 0);
    auto d = const_auto_cast(table_);
    tuple max_exp;
    for (uint64_t z = 0; z < x; ++z)
      max_exp.v[z] = u64::power(base_.v[z], window_size_);
    for (uint64_t i = 0; i < 256; ++i) {
      for (uint64_t z = 0; z < x; ++z)
        d[i][0].v[z] = u64::mod(p61 - u64::mod(i * (uint128_t)(max_exp.v[z])));
      for (uint64_t j = 1; j < 256; ++j) {
        for (uint64_t z = 0; z < x; ++z)
          d[i][j].v[z] = u64::mod(d[i][j - 1].v[z] + 1);
      }
    }
  };

  sliding_window_multi61(uint64_t const window_size)
      : sliding_window_multi61(window_size, tuple().apply([](uint64_t) {
          return u64::random(1, p61 - 1);
        })){};

  template <ByteType T>
  KRINLNFN tuple roll_right(tuple fp, T const pop_left,
                            T const push_right) const {
    auto const &lookup = table_[pop_left][push_right];
    for (uint64_t z = 0; z < x; ++z) {
      if (base_.v[z] >= p61 || fp.v[z] >= p61 || lookup.v[z] >= p61)
        __builtin_unreachable();
      else
        fp.v[z] = u64::mod(((uint128_t)base_.v[z]) * fp.v[z] + lookup.v[z]);
    }
    return fp;
  }

  template <ByteType T>
  KRINLNFN tuple roll_right(tuple fp, T const push_right) const {
    for (uint64_t z = 0; z < x; ++z) {
      if (base_.v[z] >= p61 || fp.v[z] >= p61)
        __builtin_unreachable();
      else
        fp.v[z] = u64::mod(((uint128_t)base_.v[z]) * fp.v[z] + push_right);
    }
    return fp;
  }

  inline tuple base() const { return base_; }
  inline uint64_t window_size() const { return window_size_; }
  inline uint64_t bits() const { return x * 61; }
  inline double collision_rate() const { return collision_rate_; }
};

}  // namespace u64

}  // namespace kr_fingerprinting

