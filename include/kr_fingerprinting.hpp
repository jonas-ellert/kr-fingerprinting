#pragma once

#include <bit>
#include <concepts>
#include <cstdint>
#include <ostream>
#include <random>

namespace kr_fingerprinting {

template <typename T>
concept ByteType = std::unsigned_integral<T> && (sizeof(T) == 1);

using uint128_t = unsigned __int128;

// trailing zero count
constexpr uint64_t countr_zero(uint128_t v) {
  return (v << 64) ? std::countr_zero((uint64_t)v)
                   : (std::countr_zero((uint64_t)(v >> 64)) + 64);
}

template <uint64_t s>
struct MERSENNE {
  using uintX_t = uint128_t;
  constexpr static uintX_t value = (((uintX_t)1ULL) << s) - 1;
  constexpr static uintX_t shift = s;
};

struct MERSENNE61 {
  using uintX_t = uint64_t;
  constexpr static uintX_t value = (((uintX_t)1ULL) << 61) - 1;
  constexpr static uintX_t shift = countr_zero(value + 1);
};
struct MERSENNE89 {
  using uintX_t = uint128_t;
  constexpr static uintX_t value = (((uintX_t)1ULL) << 89) - 1;
  constexpr static uintX_t shift = countr_zero(value + 1);
};
struct MERSENNE107 {
  using uintX_t = uint128_t;
  constexpr static uintX_t value = (((uintX_t)1ULL) << 107) - 1;
  constexpr static uintX_t shift = countr_zero(value + 1);
};
struct MERSENNE127 {
  using uintX_t = uint128_t;
  constexpr static uintX_t value = (((uintX_t)1ULL) << 127) - 1;
  constexpr static uintX_t shift = countr_zero(value + 1);
};

template <typename p>
concept MersennePrime =
    true || std::is_same_v<p, MERSENNE61> || std::is_same_v<p, MERSENNE89> ||
    std::is_same_v<p, MERSENNE107> || std::is_same_v<p, MERSENNE127>;

template <MersennePrime p>
struct kr_fingerprinter {
  constexpr static bool is_64_bit = std::is_same_v<p, MERSENNE61>;
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
    return (i >= prime) ? (i - prime) : i;
  }

  // this assumes a,b < prime
  inline constexpr static uint128_t mult(uintX_t const a, uintX_t const b) {
    if constexpr (is_64_bit) {
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

      uint128_t const carry =
          ((uint128_t)(uint64_t)m1 + (uint128_t)(uint64_t)m2 + (l >> 64)) >> 64;

      uint128_t const h128 = h + (m1 >> 64) + (m2 >> 64) + carry;
      uint128_t const l128 = l + (m1 << 64) + (m2 << 64);

      uint128_t sum =
          ((h128 << (128 - shift)) | (l128 >> shift)) + (l128 & prime);

      return modulo(sum);
    }
  }

  // this assumes a,b,c < prime
  // a * b + c
  inline constexpr static uint128_t mult_add(uintX_t const a, uintX_t const b,
                                             uintX_t const c) {
    if constexpr (is_64_bit) {
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

      uint128_t const carry =
          (((l >> 64) + (c >> 64) + (uint64_t)m1) + (uint64_t)m2) >> 64;

      uint128_t const h128 = h + (m1 >> 64) + (m2 >> 64) + carry;
      uint128_t const l128 = l + c + (m1 << 64) + (m2 << 64);

      uint128_t sum =
          ((h128 << (128 - shift)) | (l128 >> shift)) + (l128 & prime);

      return modulo(sum);
    }
  }

  // this assumes a < prime
  inline constexpr static uint128_t square(uintX_t const a) {
    if constexpr (is_64_bit) {
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
    if constexpr (is_64_bit)
      return modulo(mult(a, b));
    else
      return mult(a, b);
  }

  inline constexpr static uintX_t mult_add_modulo(uintX_t const a,
                                                  uintX_t const b,
                                                  uintX_t const c) {
    if constexpr (is_64_bit)
      return modulo(mult_add(a, b, c));
    else
      return mult_add(a, b, c);
  }

  inline constexpr static uintX_t square_modulo(uintX_t const a) {
    if constexpr (is_64_bit)
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

    if constexpr (is_64_bit) {
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

    inline uintX_t roll_right(uintX_t const fp, uintX_t const pop_left,
                              uintX_t const push_right) {
      constexpr static uint128_t maxproduct =
          (is_64_bit) ? (((uint128_t)prime) * prime) : ((uint128_t)prime);

      uint128_t const shifted_fingerprint = mult(base_, fp);
      uint128_t const pop = maxproduct - mult(max_exponent_, pop_left);

      if constexpr (is_64_bit) {
        return modulo(shifted_fingerprint + pop + push_right);
      } else {
        return modulo(modulo(shifted_fingerprint + pop) + push_right);
      }
    }

    // this could be more efficient!
    /*inline uintX_t roll_left(uintX_t const fp, uintX_t const push_left,
                             uintX_t const pop_right) {
      uintX_t const popped_fingerprint = prime + fp - pop_right;
      uint128_t const push = mult(max_exponent_, push_left);

      if constexpr (is_64_bit) {
        return mult_modulo(modulo(push + popped_fingerprint), inverse_base_);
      } else {
        return mult_modulo(modulo(push + modulo(popped_fingerprint)),
                           inverse_base_);
      }
    }*/

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
