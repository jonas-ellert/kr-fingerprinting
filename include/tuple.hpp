#pragma once

#include <compare>
#include <cstdint>
#include <cstring>

namespace kr_fingerprinting {

namespace kr_tuple {

template <uint64_t x>
struct tuple {
  constexpr static uint64_t size = x;
  uint64_t v[x] = {};

  template <typename T>
  tuple &apply(T const &t) {
    for (uint64_t z = 0; z < x; ++z) v[z] = t(v[z]);
    return *this;
  }

  __attribute__((always_inline)) inline bool operator==(tuple const &o) const {
    __extension__ using uint128_t = unsigned __int128;
    static_assert(x > 0);
    static_assert(sizeof(tuple) == 8 * x);
    static_assert(sizeof(uint64_t) == 8);
    static_assert(sizeof(uint128_t) == 16);

    if constexpr (x == 1) {
      return v[0] == o.v[0];

    } else if constexpr (x == 2) {
      return (*((uint128_t const *)((void const *)this))) ==
             (*((uint128_t const *)((void const *)&o)));

    } else if constexpr (x == 3) {
      return ((*((uint128_t const *)((void const *)this))) ==
              (*((uint128_t const *)((void const *)&o)))) &&
             (v[2] == o.v[2]);
    } else if constexpr (x == 4) {
      return ((*((uint128_t const *)((void const *)this))) ==
              (*((uint128_t const *)((void const *)&o)))) &&
             ((*(((uint128_t const *)((void const *)this)) + 1)) ==
              (*(((uint128_t const *)((void const *)&o)) + 1)));
    } else {
      return std::memcmp(this, &o, sizeof(tuple)) == 0;
    }
  }

  __attribute__((always_inline)) inline bool operator!=(tuple const &o) const {
    __extension__ using uint128_t = unsigned __int128;
    static_assert(x > 0);
    static_assert(sizeof(tuple) == 8 * x);
    static_assert(sizeof(uint64_t) == 8);
    static_assert(sizeof(uint128_t) == 16);

    if constexpr (x == 1) {
      return v[0] != o.v[0];

    } else if constexpr (x == 2) {
      return (*((uint128_t const *)((void const *)this))) !=
             (*((uint128_t const *)((void const *)&o)));

    } else if constexpr (x == 3) {
      return ((*((uint128_t const *)((void const *)this))) !=
              (*((uint128_t const *)((void const *)&o)))) ||
             (v[2] != o.v[2]);
    } else if constexpr (x == 4) {
      return ((*((uint128_t const *)((void const *)this))) !=
              (*((uint128_t const *)((void const *)&o)))) ||
             ((*(((uint128_t const *)((void const *)this)) + 1)) !=
              (*(((uint128_t const *)((void const *)&o)) + 1)));
    } else {
      return std::memcmp(this, &o, sizeof(tuple)) == 0;
    }
  }

  __attribute__((always_inline)) inline bool operator<(tuple const &o) const {
    __extension__ using uint128_t = unsigned __int128;
    static_assert(x > 0);
    static_assert(sizeof(tuple) == 8 * x);
    static_assert(sizeof(uint64_t) == 8);
    static_assert(sizeof(uint128_t) == 16);

    if constexpr (x == 1) {
      return v[0] < o.v[0];

    } else if constexpr (x == 2) {
      return (*((uint128_t const *)((void const *)this))) <
             (*((uint128_t const *)((void const *)&o)));

    } else if constexpr (x == 3) {
      return ((*((uint128_t const *)((void const *)this))) <
              (*((uint128_t const *)((void const *)&o)))) ||
             (((*((uint128_t const *)((void const *)this))) ==
               (*((uint128_t const *)((void const *)&o)))) &&
              (v[2] < o.v[2]));

    } else if constexpr (x == 4) {
      return ((*((uint128_t const *)((void const *)this))) <
              (*((uint128_t const *)((void const *)&o)))) ||
             (((*((uint128_t const *)((void const *)this))) ==
               (*((uint128_t const *)((void const *)&o)))) &&
              ((*(((uint128_t const *)((void const *)this)) + 1)) <
               (*(((uint128_t const *)((void const *)&o)) + 1))));

    } else {
      return std::memcmp(this, &o, sizeof(tuple)) < 0;
    }
  }

  __attribute__((always_inline)) inline bool operator<=(tuple const &o) const {
    __extension__ using uint128_t = unsigned __int128;
    static_assert(x > 0);
    static_assert(sizeof(tuple) == 8 * x);
    static_assert(sizeof(uint64_t) == 8);
    static_assert(sizeof(uint128_t) == 16);

    if constexpr (x == 1) {
      return v[0] <= o.v[0];

    } else if constexpr (x == 2) {
      return (*((uint128_t const *)((void const *)this))) <=
             (*((uint128_t const *)((void const *)&o)));

    } else if constexpr (x == 3) {
      return ((*((uint128_t const *)((void const *)this))) <
              (*((uint128_t const *)((void const *)&o)))) ||
             (((*((uint128_t const *)((void const *)this))) ==
               (*((uint128_t const *)((void const *)&o)))) &&
              (v[2] <= o.v[2]));

    } else if constexpr (x == 4) {
      return ((*((uint128_t const *)((void const *)this))) <
              (*((uint128_t const *)((void const *)&o)))) ||
             (((*((uint128_t const *)((void const *)this))) ==
               (*((uint128_t const *)((void const *)&o)))) &&
              ((*(((uint128_t const *)((void const *)this)) + 1)) <=
               (*(((uint128_t const *)((void const *)&o)) + 1))));

    } else {
      return std::memcmp(this, &o, sizeof(tuple)) <= 0;
    }
  }

  bool operator>(tuple const &o) const { return o < *this; }
  bool operator>=(tuple const &o) const { return o <= *this; }
};

}  // namespace kr_tuple

}  // namespace kr_fingerprinting

namespace std {

template <uint64_t x>
std::string to_string(kr_fingerprinting::kr_tuple::tuple<x> const &t) {
  std::stringstream str;
  str << t.v[0];
  for (uint64_t i = 1; i < x; ++i) str << "-" << t.v[i];
  return str.str();
}

template <uint64_t x>
std::ostream &operator<<(std::ostream &s,
                         kr_fingerprinting::kr_tuple::tuple<x> const &t) {
  s << to_string(t);
  return s;
}

}  // namespace std
