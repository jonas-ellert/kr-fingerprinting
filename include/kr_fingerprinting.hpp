#pragma once

#include <bit>
#include <concepts>
#include <cstdint>
#include <ostream>
#include <random>

namespace kr_fingerprinting {

template <typename T>
concept ByteType = std::unsigned_integral<T> && (sizeof(T) == 1);

__extension__ using uint128_t = unsigned __int128;

// trailing zero count
constexpr uint64_t countr_zero(uint128_t v) {
  return (v << 64) ? std::countr_zero((uint64_t)v)
                   : (std::countr_zero((uint64_t)(v >> 64)) + 64);
}

template <uint64_t s>
struct TWO_POW_MINUS_ONE {
  using uintX_t = std::conditional_t<(s > 63), uint128_t, uint64_t>;
  constexpr static uintX_t value = (((uintX_t)1ULL) << s) - 1;
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

  // this assumes a,b < prime
  inline constexpr static uint128_t mult(uintX_t const a, uintX_t const b) {
    if constexpr (shift < 64) {
      return ((uint128_t)a) * b;
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
        // this only works because we have sufficiently many overflow bits
        uint128_t const m = m1 + m2;
        uint128_t const sum = (l & prime) + (l >> shift) + ((m << 64) & prime) +
                              (m >> (shift - 64)) +
                              ((h << (128 - shift)) & prime) +
                              (h >> (2 * shift - 128));
        return modulo(sum);
      } else {
        uint128_t const carry =
            (((l >> 64) + (uint64_t)m1) + (uint64_t)m2) >> 64;

        uint128_t const h128 = h + (m1 >> 64) + (m2 >> 64) + carry;
        uint128_t const l128 = l + (m1 << 64) + (m2 << 64);

        uint128_t sum =
            ((h128 << (128 - shift)) | (l128 >> shift)) + (l128 & prime);

        return modulo(sum);
      }
    }
  }

  // this assumes a,b,c < prime
  // a * b + c
  inline constexpr static uint128_t mult_add(uintX_t const a, uintX_t const b,
                                             uintX_t const c) {
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

  // this assumes a < prime
  inline constexpr static uint128_t square(uintX_t const a) {
    if constexpr (shift < 64) {
      return ((uint128_t)a) * a;
    } else {
      uint128_t const al = (uint64_t)a;
      uint128_t const ah = a >> 64;

      uint128_t const h = ah * ah;
      uint128_t const m = ah * al;
      uint128_t const l = al * al;

      uint128_t const carry = (((m << 64) >> 63) + (l >> 64)) >> 64;

      uint128_t const h128 = h + ((m >> 64) << 1) + carry;
      uint128_t const l128 = l + (m << 65);

      uint128_t sum =
          ((h128 << (128 - shift)) | (l128 >> shift)) + (l128 & prime);

      return modulo(sum);
    }
  }

  inline constexpr static uintX_t mult_modulo(uintX_t const a,
                                              uintX_t const b) {
    if constexpr (shift < 64)
      return modulo(mult(a, b));
    else
      return mult(a, b);
  }

  inline constexpr static uintX_t mult_add_modulo(uintX_t const a,
                                                  uintX_t const b,
                                                  uintX_t const c) {
    if constexpr (shift < 64)
      return modulo(mult_add(a, b, c));
    else
      return mult_add(a, b, c);
  }

  inline constexpr static uintX_t square_modulo(uintX_t const a) {
    if constexpr (shift < 64)
      return modulo(square(a));
    else
      return square(a);
  }

  // fast squaring
  inline constexpr static uintX_t power(uintX_t base, uintX_t exponent) {
    uintX_t result = 1;
    while (exponent > 0) {
      if (exponent & 1ULL) result = mult_modulo(base, result);
      base = square_modulo(base);
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
      uintX_t const result = d_lo(g) | (uint128_t)(d_hi(g));
      return (!result) ? random_base() : result;
    }
  }

  inline static uintX_t inverse_base(uintX_t const base) {
    // using eulers theorem
    return power(base, prime - 2);
  }

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

  template <bool large_precomputation>
  struct sliding_window_precompute {
   private:
    uint64_t const window_size_;
    uintX_t const base_;
    uintX_t const inverse_base_ = power(base_, prime - 2);
    uintX_t const max_exponent_ = power(base_, window_size_);

    using lookup =
        std::conditional_t<large_precomputation, push_pop_lookup, pop_lookup>;
    lookup const table_ = lookup(max_exponent_);

   public:
    sliding_window_precompute(uint64_t const window_size, uintX_t const base)
        : window_size_(window_size), base_(modulo(base)){};

    template <ByteType T>
    inline uintX_t roll_right(uintX_t const fp, T const pop_left,
                              T const push_right) {
      if constexpr (large_precomputation) {
        return mult_add_modulo(base_, fp, table_(pop_left, push_right));
      } else {
        return modulo(mult_add(base_, fp, table_(pop_left)) + push_right);
      }
    }

    template <ByteType T>
    inline uintX_t roll_right(uintX_t const fp, T const push_right) {
      return mult_add_modulo(base_, fp, push_right);
    }

    inline uintX_t base() const { return base_; }
  };

  using sliding_window = sliding_window_precompute<true>;
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
