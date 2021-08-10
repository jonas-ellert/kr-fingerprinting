#pragma once

#include <bit>
#include <limits>
#include <random>
#include <type_traits>

namespace kr_fingerprinting {

__extension__ using uint128_t = unsigned __int128;

template <typename T>
concept ByteType = std::unsigned_integral<T> && (sizeof(T) == 1);

template <typename T1, typename T2>
auto &cast(T2 &t) {
  static_assert(sizeof(T1) == sizeof(T2));
  return *reinterpret_cast<T1 *>(&t);
}

template <typename T>
struct table2d {
 private:
  T data[256][256];

 public:
  template <typename initializer>
  table2d(initializer &&init) {
    init(data);
  }

  template <ByteType B>
  inline T operator()(B const l, B const r) const {
    return data[l][r];
  }
};

namespace u64 {

template <uint64_t modulus>
inline constexpr uint64_t mod(uint128_t const value) {
  static_assert(std::popcount(modulus) < 64);
  static_assert(std::popcount(modulus) == std::countr_one(modulus));
  constexpr uint64_t shift = std::popcount(modulus);
  uint64_t const i = (value & modulus) + (value >> shift);
  return (i & modulus) + (i >> shift);
}

// fast squaring
template <uint64_t modulus>
inline constexpr uint64_t power(uint64_t base, uint64_t exponent) {
  uint64_t result = 1;
  uint128_t b = base;
  while (exponent > 0) {
    if (exponent & 1ULL) result = mod<modulus>(b * result);
    b = mod<modulus>(b * b);
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
  // mersenne prime 2^61 - 1
  constexpr static uint64_t s = 61;
  constexpr static uint64_t p = (1ULL << s) - 1;

  uint64_t const window_size_;
  uint64_t const base_;

  double const collision_rate_ = ((double)window_size_ - 1) / p;

  table2d<uint64_t> const table_ = table2d<uint64_t>([&](uint64_t d[][256]) {
    uint64_t const max_exponent = u64::power<p>(base_, window_size_);
    for (uint64_t i = 0; i < 256; ++i) {
      d[i][0] = u64::mod<p>(p - u64::mod<p>(i * (uint128_t)max_exponent));
      for (uint64_t j = 1; j < 256; ++j) {
        d[i][j] = u64::mod<p>(d[i][j - 1] + 1);
      }
    }
  });

 public:
  using fingerprint_type = uint64_t;

  sliding_window61(uint64_t const window_size, uint64_t const base)
      : window_size_(window_size), base_(u64::mod<p>(base)){};

  sliding_window61(uint64_t const window_size)
      : sliding_window61(window_size, u64::random(1, p - 1)){};

  template <ByteType T>
  inline uint64_t roll_right(uint64_t const fp, T const pop_left,
                             T const push_right) const {
    auto lookup = table_(pop_left, push_right);
    if (base_ >= p || fp >= p || lookup >= p)
      __builtin_unreachable();
    else
      return u64::mod<p>(((uint128_t)base_) * fp + lookup);
  }

  template <ByteType T>
  inline uint64_t roll_right(uint64_t const fp, T const push_right) const {
    if (base_ >= p || fp >= p)
      __builtin_unreachable();
    else
      return u64::mod<p>(((uint128_t)base_) * fp + push_right);
  }

  inline uint64_t base() const { return base_; }
  inline uint64_t window_size() const { return window_size_; }
  inline uint64_t bits() const { return s; }
  inline double collision_rate() const { return collision_rate_; }
};

template <uint64_t x>
struct sliding_window_multi61 {
 public:
  struct tuple {
    uint64_t v[x] = {};
    inline uint64_t &operator[](uint64_t i) { return v[i]; }
    inline uint64_t const &operator[](uint64_t i) const { return v[i]; }
    tuple &mod() {
      for (uint64_t z = 0; z < x; ++z) v[z] = u64::mod<p>(v[z]);
      return *this;
    }
    tuple &rnd() {
      for (uint64_t z = 0; z < x; ++z) v[z] = u64::random(1, p - 1);
      return *this;
    }
    auto operator<=>(tuple const &) const = default;
  };

 private:
  // mersenne prime 2^61 - 1
  constexpr static uint64_t s = 61;
  constexpr static uint64_t p = (1ULL << s) - 1;

  uint64_t const window_size_;
  tuple const base_;
  double const collision_rate_ = std::pow(((double)window_size_ - 1) / p, x);

  table2d<tuple> const table_ = table2d<tuple>([&](tuple d[][256]) {
    tuple max_exp;
    for (uint64_t z = 0; z < x; ++z)
      max_exp[z] = u64::power<p>(base_[z], window_size_);
    for (uint64_t i = 0; i < 256; ++i) {
      for (uint64_t z = 0; z < x; ++z)
        d[i][0][z] = u64::mod<p>(p - u64::mod<p>(i * (uint128_t)(max_exp[z])));
      for (uint64_t j = 1; j < 256; ++j) {
        for (uint64_t z = 0; z < x; ++z)
          d[i][j][z] = u64::mod<p>(d[i][j - 1][z] + 1);
      }
    }
  });

  constexpr static bool is_uint128_t = sizeof(tuple) == sizeof(uint128_t);

 public:
  using fingerprint_type = std::conditional_t<is_uint128_t, uint128_t, tuple>;

  sliding_window_multi61(uint64_t const window_size, fingerprint_type base)
      : window_size_(window_size), base_((cast<tuple>(base)).mod()) {
    static_assert(x >= 2);
  };

  sliding_window_multi61(uint64_t const window_size)
      : sliding_window_multi61(window_size,
                               cast<fingerprint_type>(tuple().rnd())){};

  template <ByteType T>
  inline fingerprint_type roll_right(fingerprint_type fingerprint,
                                     T const pop_left,
                                     T const push_right) const {
    tuple &fp = cast<tuple>(fingerprint);
    auto const &lookup = table_(pop_left, push_right);
    for (uint64_t z = 0; z < x; ++z) {
      if (base_[z] >= p || fp[z] >= p || lookup[z] >= p)
        __builtin_unreachable();
      else
        fp[z] = u64::mod<p>(((uint128_t)base_[z]) * fp[z] + lookup[z]);
    }
    return fingerprint;
  }

  template <ByteType T>
  inline fingerprint_type roll_right(fingerprint_type fingerprint,
                                     T const push_right) const {
    tuple &fp = cast<tuple>(fingerprint);
    for (uint64_t z = 0; z < x; ++z) {
      if (base_[z] >= p || fp[z] >= p)
        __builtin_unreachable();
      else
        fp[z] = u64::mod<p>(((uint128_t)base_[z]) * fp[z] + push_right);
    }
    return fingerprint;
  }

  inline fingerprint_type base() const {
    return cast<fingerprint_type const>(base_);
  }
  inline uint64_t window_size() const { return window_size_; }
  inline uint64_t bits() const { return s * x; }
  inline double collision_rate() const { return collision_rate_; }
};

}  // namespace u64

}  // namespace kr_fingerprinting
