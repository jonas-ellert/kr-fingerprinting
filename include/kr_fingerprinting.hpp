#pragma once

#include <bit>
#include <concepts>
#include <cstdint>
#include <ostream>
#include <random>

namespace kr_fingerprinting {

__extension__ using uint128_t = unsigned __int128;

template <typename T>
concept ByteType = std::unsigned_integral<T> && (sizeof(T) == 1);

// trailing zero count
constexpr uint64_t countr_zero(uint128_t v) {
  return (v << 64) ? std::countr_zero((uint64_t)v)
                   : (std::countr_zero((uint64_t)(v >> 64)) + 64);
}

template <uint64_t s>
struct TWO_POW_MINUS_ONE {
  using uintX_t = std::conditional_t<(s > 64), uint128_t, uint64_t>;
  constexpr static uintX_t value = (((uint128_t)1) << s) - 1;
  constexpr static uintX_t shift = s;
};

constexpr bool is_mersenne_power(uint64_t const p) {
  constexpr uint64_t powers[] = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127};
  for (uint64_t i = 0; i < sizeof(powers) / sizeof(uint64_t); ++i)
    if (powers[i] == p) return true;
  return false;
}

using MERSENNE61 = TWO_POW_MINUS_ONE<61>;
using MERSENNE89 = TWO_POW_MINUS_ONE<89>;
using MERSENNE107 = TWO_POW_MINUS_ONE<107>;
using MERSENNE127 = TWO_POW_MINUS_ONE<127>;

template <typename T>
concept MersennePrime = is_mersenne_power(T::shift) &&
                        std::is_same_v<T, TWO_POW_MINUS_ONE<T::shift>>;

template <MersennePrime p>
struct kr_fingerprinter {
  using uintX_t = p::uintX_t;
  constexpr static uintX_t prime = p::value;
  constexpr static uintX_t shift = p::shift;

  inline constexpr static uintX_t modulo(uint128_t const value) {
    /*if constexpr (is_64_bit) {
      // this assumes value < 2^122
      uint128_t const v = value + 1;
      uint64_t const z = ((v >> shift) + v) >> shift;
      return (value + z) & prime;
    } else {
      // this assumes value < 2^128 - 1
      uint128_t const z = (value + 1) >> shift;
      return (value + z) & prime;
    }*/

    // another option for any type:
    // this assumes value < prime^2
    uintX_t const i = (value & prime) + (value >> shift);
    if constexpr (shift < 64)
      return (i & prime) + (i >> shift);
    else
      return (i >= prime) ? (i - prime) : i;
  }

  // this assumes a,b,c < prime
  // a * b + c
  inline constexpr static uint128_t mult_add(uintX_t const a, uintX_t const b,
                                             uintX_t const c) {
    if (a >= prime) __builtin_unreachable();
    if (b >= prime) __builtin_unreachable();
    if (c >= prime) __builtin_unreachable();
    if constexpr (shift < 64) {
      return (((uint128_t)a) * b) + c;
    } else {
      uint128_t const al = (uint64_t)a;
      uint128_t const ah = a >> 64;
      uint128_t const bl = (uint64_t)b;
      uint128_t const bh = b >> 64;

      uint128_t const h = ah * bh;
      uint128_t const m1 = ah * bl;
      uint128_t const m2 = bh * al;
      uint128_t const l = al * bl;

      if constexpr (shift < 127) {
        if (h >= prime) __builtin_unreachable();
        if (l >= prime) __builtin_unreachable();
        if (m1 >= prime) __builtin_unreachable();
        if (m2 >= prime) __builtin_unreachable();
        // this only works because we have sufficiently many overflow bits
        uint128_t const m = m1 + m2;
        uint128_t const sum = c + (l & prime) + (l >> shift) +
                              ((m << 64) & prime) + (m >> (shift - 64)) +
                              ((h << (128 - shift)) & prime) +
                              (h >> (2 * shift - 128));
        return modulo(sum);
      } else {
        // for p127 the single overflow bit is not sufficient
        uint128_t const carry =
            ((l >> 64) + (c >> 64) + (uint64_t)m1 + (uint64_t)m2) >> 64;

        uint128_t const h128 = h + (m1 >> 64) + (m2 >> 64) + carry;
        uint128_t const l128 = l + c + (m1 << 64) + (m2 << 64);

        uint128_t sum =
            ((h128 << (128 - shift)) | (l128 >> shift)) + (l128 & prime);

        return modulo(sum);
      }
    }
  }

  inline constexpr static uint128_t mult(uintX_t const a, uintX_t const b) {
    return mult_add(a, b, 0);
  }

  inline constexpr static uintX_t mult_add_modulo(uintX_t const a,
                                                  uintX_t const b,
                                                  uintX_t const c) {
    if constexpr (shift < 64)
      return modulo(mult_add(a, b, c));
    else
      return mult_add(a, b, c);
  }

  inline constexpr static uintX_t mult_modulo(uintX_t const a,
                                              uintX_t const b) {
    return mult_add_modulo(a, b, 0);
  }

  // fast squaring
  inline constexpr static uintX_t power(uintX_t base, uintX_t exponent) {
    uintX_t result = 1;
    while (exponent > 0) {
      if (exponent & 1ULL) result = mult_modulo(base, result);
      base = mult_modulo(base, base);
      exponent >>= 1;
    }
    return result;
  }

  inline static uintX_t random_base() {
    using unif_dist = std::uniform_int_distribution<uint64_t>;
    static std::random_device seed;
    static std::mt19937_64 g = std::mt19937_64(seed());

    if constexpr (shift < 64) {
      static unif_dist d(1, prime - 2);
      return d(g);
    } else {
      static unif_dist d_hi(0, (uint64_t)((prime - 2) >> 64));
      static unif_dist d_lo(std::numeric_limits<uint64_t>::min(),
                            std::numeric_limits<uint64_t>::max());
      // avoid returning 0
      uintX_t const result = d_lo(g) | (((uint128_t)d_hi(g)) << 64);
      return (!result) ? random_base() : result;
    }
  }

  inline static uintX_t inverse_base(uintX_t const base) {
    // using eulers theorem
    return power(base, prime - 2);
  }

  template <bool large_precomputation>
  struct sliding_window_precompute {
   private:
    struct pop_lookup {
     private:
      uintX_t data[256];

     public:
      pop_lookup(uintX_t const base) {
        for (uint64_t i = 0; i < 256; ++i) {
          data[i] = prime - mult_modulo(i, base);
        }
      }

      template <ByteType T>
      inline uintX_t operator()(T const t) const {
        return data[t];
      }
    };

    struct push_pop_lookup {
     private:
      uintX_t data[256][256];

     public:
      push_pop_lookup(uintX_t const base) {
        for (uint64_t i = 0; i < 256; ++i) {
          data[i][0] = prime - mult_modulo(i, base);
          for (uint64_t j = 1; j < 256; ++j) {
            data[i][j] = modulo(data[i][j - 1] + 1);
          }
        }
      }

      template <ByteType T>
      inline uintX_t operator()(T const l, T const r) const {
        return data[l][r];
      }
    };

    uint64_t const window_size_;
    uintX_t const base_;
    uintX_t const max_exponent_ = power(base_, window_size_);

    using lookup =
        std::conditional_t<large_precomputation, push_pop_lookup, pop_lookup>;
    lookup const table_ = lookup(max_exponent_);

   public:
    sliding_window_precompute(uint64_t const window_size, uintX_t const base)
        : window_size_(window_size), base_(modulo(base)){};

    template <ByteType T>
    inline uintX_t roll_right(uintX_t const fp, T const pop_left,
                              T const push_right) const {
      if constexpr (large_precomputation) {
        return mult_add_modulo(base_, fp, table_(pop_left, push_right));
      } else {
        return modulo(mult_add(base_, fp, table_(pop_left)) + push_right);
      }
    }

    template <ByteType T>
    inline uintX_t roll_right(uintX_t const fp, T const push_right) const {
      return mult_add_modulo(base_, fp, push_right);
    }

    inline uintX_t base() const { return base_; }
    inline uint64_t window_size() const { return window_size_; }
    inline uint64_t bits() const { return shift; }
  };

  using sliding_window = sliding_window_precompute<true>;
};

using sliding_window61 = kr_fingerprinter<MERSENNE61>::sliding_window;
using sliding_window89 = kr_fingerprinter<MERSENNE89>::sliding_window;
using sliding_window107 = kr_fingerprinter<MERSENNE107>::sliding_window;
using sliding_window127 = kr_fingerprinter<MERSENNE127>::sliding_window;

inline static uint128_t random_base_pair61() {
  return (((uint128_t)kr_fingerprinter<MERSENNE61>::random_base()) << 64) |
         kr_fingerprinter<MERSENNE61>::random_base();
}

struct sliding_window122 {
 private:
  using F = kr_fingerprinter<MERSENNE61>;

  struct pair {
    uint64_t a;
    uint64_t b;
  };

  struct push_pop_lookup {
   private:
    pair data[256][256];

   public:
    push_pop_lookup(uint64_t const max_exp1, uint64_t const max_exp2) {
      for (uint64_t i = 0; i < 256; ++i) {
        data[i][0].a = F::modulo(F::prime - F::mult_modulo(i, max_exp1));
        data[i][0].b = F::modulo(F::prime - F::mult_modulo(i, max_exp2));
        for (uint64_t j = 1; j < 256; ++j) {
          data[i][j].a = F::modulo(data[i][j - 1].a + 1);
          data[i][j].b = F::modulo(data[i][j - 1].b + 1);
        }
      }
    }

    template <ByteType T>
    inline pair operator()(T const l, T const r) const {
      return data[l][r];
    }
  };

  uint64_t const window_size_;
  uint64_t const base1_;
  uint64_t const base2_;
  uint64_t const max_exponent1_ = F::power(base1_, window_size_);
  uint64_t const max_exponent2_ = F::power(base2_, window_size_);
  push_pop_lookup const table_ =
      push_pop_lookup(max_exponent1_, max_exponent2_);

 public:
  sliding_window122(uint64_t const window_size, uint64_t const base1,
                    uint64_t const base2)
      : window_size_(window_size),
        base1_(F::modulo(base1)),
        base2_(F::modulo(base2)){};

  sliding_window122(uint64_t const window_size, uint128_t const basepair)
      : sliding_window122(window_size, basepair >> 64, (uint64_t)basepair){};

  template <ByteType T>
  inline uint128_t roll_right(uint128_t const fp, T const pop_left,
                              T const push_right) const {
    auto const lookup = table_(pop_left, push_right);
    uint128_t const l = F::mult_add_modulo(base1_, fp >> 64, lookup.a);
    uint64_t const r = F::mult_add_modulo(base2_, (uint64_t)fp, lookup.b);
    return (l << 64) | r;
  }

  template <ByteType T>
  inline uint128_t roll_right(uint128_t const fp, T const push_right) const {
    return roll_right(fp, (T)0, push_right);
  }

  inline uint64_t base1() const { return base1_; }
  inline uint64_t base2() const { return base2_; }
  inline uint128_t base() const { return ((uint128_t)base1_) << 64 | base2_; }
  inline uint64_t window_size() const { return window_size_; }
  inline uint64_t bits() const { return 122; }
};

}  // namespace kr_fingerprinting

inline std::ostream &operator<<(std::ostream &out,
                                kr_fingerprinting::uint128_t value) {
  std::ostream::sentry s(out);
  if (s) {
    char buffer[64];  // 39 should be enough
    char *digit = &(buffer[64]);
    do {
      *(--digit) = "0123456789"[value % 10];
      value /= 10;
    } while (value != 0);
    int len = &(buffer[64]) - digit;
    if (out.rdbuf()->sputn(digit, len) != len) {
      out.setstate(std::ios_base::badbit);
    }
  }
  return out;
}

namespace std {
std::string to_string(kr_fingerprinting::uint128_t value) {
  char buffer[64];  // 39 should be enough
  char *digit = &(buffer[64]);
  do {
    *(--digit) = "0123456789"[value % 10];
    value /= 10;
  } while (value != 0);
  return std::string(digit, &buffer[64]);
}
}  // namespace std
