#pragma once

#include "kr-fingerprinting128.hpp"

namespace kr_fingerprinting {
    
using sliding_window89 = u128::sliding_windowX<89>;
using sliding_window107 = u128::sliding_windowX<107>;
using sliding_window127 = u128::sliding_windowX<127>;
    
using sliding_window61 = u64::sliding_window61;
using sliding_window122 = u64::sliding_window122;
using sliding_window122b = u64::sliding_window_multi61<2>;
using sliding_window183 = u64::sliding_window_multi61<3>;
using sliding_window244 = u64::sliding_window_multi61<4>;

}  // namespace kr_fingerprinting


