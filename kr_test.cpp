#include <chrono>
#include <concepts>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

//#define inline __attribute__((always_inline)) inline

#include "include/kr-fingerprinting.hpp"
#include "include/kr-fingerprinting128.hpp"

#include "../rk-fingerprint/rolling_hash/rk_prime.hpp"

using namespace kr_fingerprinting;

struct {
  std::chrono::steady_clock::time_point begin;

  void start() { begin = std::chrono::steady_clock::now(); }

  uint64_t stop() const {
    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
        .count();
  }

  double mibs(uint64_t const ms, uint64_t const n) const {
    double const s = ms / 1000.0;
    double const mib = n / 1024.0 / 1024.0;
    return mib / s;
  }

} timer;

template <typename window_type>
KRINLNFN void mainp(std::vector<uint8_t> const &string, window_type const &w) {
  uint64_t const n = string.size();
  uint64_t const tau = w.window_size();

  // auto b = w.base();
  // std::cout << "b=" << b << std::endl;
  // std::cout << "b=" << (b >> 64) << "," << (uint64_t)b << std::endl;

  auto col = [](double rate) {
    std::stringstream stream;
    stream << "1/2^" << std::setprecision(2) << std::fixed
           << (std::floor(100 * std::log2(1.0 / rate)) / 100.0);
    return stream.str();
  };

  using uintX_t = window_type::fingerprint_type;

  std::string s = std::string("FP-") + std::to_string(w.bits());
  std::cout << s << " start!" << std::endl;
  std::cout << s << " collision rate: " << col(w.collision_rate()) << std::endl;
  timer.start();
  auto fp = uintX_t();
  for (size_t i = 0; i < tau; i++) {
    fp = w.roll_right(fp, string[i]);
  }
  for (size_t i = 0; i < n - tau; i++) {
    fp = w.roll_right(fp, string[i], string[i + tau]);
  }
  auto time = timer.stop();
  std::cout << s << " time: " << time << "[ms]"
            << " = " << timer.mibs(time, n) << "mibs" << std::endl;
  if constexpr (std::is_same_v<uintX_t, kr_fingerprinting::uint128_t>) {
    std::cout << s << " hash: " << (uint64_t)(fp >> 64) << ", " << (uint64_t)fp
              << std::endl;
  } else {
    std::cout << s << " hash: " << fp << std::endl;
  }

  uintX_t fptest = uintX_t();
  for (size_t i = n - tau; i < n; i++) {
    fptest = w.roll_right(fptest, string[i]);
  }
  std::cout << s << " correct=" << (fptest == fp) << std::endl;
}

template <typename window_type>
KRINLNFN void mainp1(std::vector<uint8_t> const &string, window_type const &w) {
  uint64_t const n = string.size();
  uint64_t const tau = w.window_size();

  // auto b = w.base();
  // std::cout << "b=" << b << std::endl;
  // std::cout << "b=" << (b >> 64) << "," << (uint64_t)b << std::endl;

  auto col = [](double rate) {
    std::stringstream stream;
    stream << "1/2^" << std::setprecision(2) << std::fixed
           << (std::floor(100 * std::log2(1.0 / rate)) / 100.0);
    return stream.str();
  };

  using uintX_t = decltype(w.base());

  std::string s = std::string("FP-NOINLINE-") + std::to_string(w.bits());
  std::cout << s << " start!" << std::endl;
  std::cout << s << " collision rate: " << col(w.collision_rate()) << std::endl;
  timer.start();
  auto fp = uintX_t();
  for (size_t i = 0; i < tau; i++) {
    fp = w.roll_right(fp, string[i]);
  }
  for (size_t i = 0; i < n - tau; i++) {
    fp = kr_fingerprinting::roll(w, fp, string[i], string[i + tau]);
  }
  auto time = timer.stop();
  std::cout << s << " time: " << time << "[ms]"
            << " = " << timer.mibs(time, n) << "mibs" << std::endl;
  // std::cout << s << " hash: " << result << std::endl;

  uintX_t fptest = uintX_t();
  for (size_t i = n - tau; i < n; i++) {
    fptest = w.roll_right(fptest, string[i]);
  }
  std::cout << s << " correct=" << (fptest == fp) << std::endl;
}

template <uint64_t shift>
KRINLNFN void mainp2(std::vector<uint8_t> const &string, uint64_t const tau) {
  auto base = kr_fingerprinting::u64::random(0, (1ULL << 19) - 1);
  auto rk = herlez::rolling_hash::rk_prime<decltype(string.cbegin()), shift>(
      string.cbegin(), tau, base);

  size_t last_window_index = string.size() - tau;

  std::string s = std::string("ALEX") + std::to_string(shift);
  std::cout << s << " start!" << std::endl;
  timer.start();
  auto fp = rk.get_currect_fp();
  for (size_t i = 0; i < last_window_index; ++i) {
    fp = rk.roll();
  }
  auto time = timer.stop();
  std::cout << s << " time: " << time << "[ms]"
            << " = " << timer.mibs(time, string.size()) << "mibs" << std::endl;
  // std::cout << s << " hash: " << result << std::endl;

  auto rk_test =
      herlez::rolling_hash::rk_prime<decltype(string.cbegin()), shift>(
          string.cbegin() + last_window_index, tau, base);
  auto fp_test = rk_test.get_currect_fp();
  std::cout << s << " correct=" << (fp_test == fp) << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc < 2) return -1;
  uint64_t const tau = std::stoi(argv[1]);

  std::vector<uint8_t> string;

  if (argc > 2) {
    std::ifstream t(argv[2]);

    t.seekg(0, std::ios::end);
    string = std::vector<uint8_t>();
    string.reserve(t.tellg());
    t.seekg(0, std::ios::beg);

    string.assign((std::istreambuf_iterator<char>(t)),
                  std::istreambuf_iterator<char>());
    std::cout << "String loaded: " << string.size() << std::endl;
  } else {
    static std::random_device seed;
    static std::mt19937_64 g(10);
    static std::uniform_int_distribution<uint8_t> d(1, 255);
    uint64_t const n = 1ULL * 128 * 1024 * 1024;
    string.resize(n);
    for (size_t i = 0; i < n; ++i) {
      string[i] = d(g);
    }
    std::cout << "String generated." << std::endl;
  }

  mainp(string, kr_fingerprinting::sliding_window<61>(tau));
  mainp(string, kr_fingerprinting::sliding_window<122>(tau));
  mainp(string, kr_fingerprinting::sliding_window<183>(tau));
  mainp(string, kr_fingerprinting::sliding_window<244>(tau));

  mainp(string, kr_fingerprinting::sliding_window<89>(tau));
  mainp(string, kr_fingerprinting::sliding_window<107>(tau));
  mainp(string, kr_fingerprinting::sliding_window<127>(tau));

  mainp1(string, kr_fingerprinting::sliding_window<61>(tau));
  mainp1(string, kr_fingerprinting::sliding_window<122>(tau));
  mainp1(string, kr_fingerprinting::sliding_window<183>(tau));
  mainp1(string, kr_fingerprinting::sliding_window<244>(tau));

  mainp1(string, kr_fingerprinting::sliding_window<89>(tau));
  mainp1(string, kr_fingerprinting::sliding_window<107>(tau));
  mainp1(string, kr_fingerprinting::sliding_window<127>(tau));

  mainp2<61>(string, tau);
  mainp2<89>(string, tau);
  mainp2<107>(string, tau);

  return 0;
}
