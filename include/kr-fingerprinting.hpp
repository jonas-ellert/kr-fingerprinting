#pragma once

#include "kr-fingerprinting128.hpp"

namespace kr_fingerprinting {

using sliding_window61 = u64::sliding_window61;
using sliding_window122 = u64::sliding_window_multi61<2>;
using sliding_window183 = u64::sliding_window_multi61<3>;
using sliding_window244 = u64::sliding_window_multi61<4>;

using sliding_window89 = u128::sliding_windowX<89>;
using sliding_window107 = u128::sliding_windowX<107>;
using sliding_window127 = u128::sliding_windowX<127>;

template <uint64_t shift>
constexpr static auto shift_to_pointer_type() {
  if constexpr (shift == 61)
    return (sliding_window61 *)0;
  else if constexpr (shift == 122)
    return (sliding_window122 *)0;
  else if constexpr (shift == 183)
    return (sliding_window183 *)0;
  else if constexpr (shift == 244)
    return (sliding_window244 *)0;

  else if constexpr (shift == 89)
    return (sliding_window89 *)0;
  else if constexpr (shift == 107)
    return (sliding_window107 *)0;
  else if constexpr (shift == 127)
    return (sliding_window127 *)0;

  else {
    constexpr bool invalid_shift = false && shift;  // always false
    static_assert(
        invalid_shift,
        "\n\n"
        "******************************************************************\n"
        "Supported parameters for kr_fingerprinting::sliding_window<shift>:\n"
        "61,89,107,122,127,183,244\n"
        "******************************************************************\n");
    return (void *)0;
  };
}

template <uint64_t shift>
using shift_to_type =
    std::remove_pointer_t<decltype(shift_to_pointer_type<shift>())>;

template <uint64_t shift>
constexpr static bool shift_supported =
    !std::is_same_v<shift_to_type<shift>, void>;

template <uint64_t shift>
using sliding_window =
    std::enable_if_t<shift_supported<shift>, shift_to_type<shift>>;

}  // namespace kr_fingerprinting
