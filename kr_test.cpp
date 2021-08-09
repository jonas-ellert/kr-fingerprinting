#include <chrono>
#include <concepts>
#include <fstream>
#include <iostream>
#include <random>

#include "include/kr_fingerprinting.hpp"

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
void mainp(std::vector<uint8_t> const &string, window_type const &w) {
  uint64_t const n = string.size();
  uint64_t const tau = w.window_size();
  auto const b = w.base();
  std::cout << "b=" << b << std::endl;
  std::cout << "b=" << (b >> 64) << "," << (uint64_t)b << std::endl;

  auto measure = [&]<typename T>(std::string s, T && t) {
    std::cout << s << " start!" << std::endl;
    timer.start();
    auto result = t();
    auto time = timer.stop();
    std::cout << s << " time: " << time << "[ms]"
              << " = " << timer.mibs(time, n) << "mibs" << std::endl;
    std::cout << s << " hash: " << result << std::endl;
    return result;
  };

  using uintX_t = decltype(w.roll_right(0, (uint8_t)0));

  std::string name = std::string("FP") + std::to_string(w.bits());
  auto fp = measure(name, [&]() {
    uintX_t fp = 0;
    for (size_t i = 0; i < tau; i++) {
      fp = w.roll_right(fp, string[i]);
    }
    for (size_t i = 0; i < n - tau; i++) {
      fp = w.roll_right(fp, string[i], string[i + tau]);
    }
    return fp;
  });

  uintX_t fptest = 0;
  for (size_t i = n - tau; i < n; i++) {
    fptest = w.roll_right(fptest, string[i]);
  }
  std::cout << name << " correct=" << (fptest == fp) << std::endl;
}

template <uint64_t s>
using fp_type = kr_fingerprinting::kr_fingerprinter<TWO_POW_MINUS_ONE<s>>;

template <uint64_t s>
using sw_type = fp_type<s>::sliding_window;

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

  mainp(string, sw_type<61>(tau, fp_type<61>::random_base()));
  mainp(string, kr_fingerprinting::sliding_window122(
                    tau, kr_fingerprinting::random_base_pair61()));
  mainp(string, sw_type<89>(tau, fp_type<89>::random_base()));
  mainp(string, sw_type<107>(tau, fp_type<107>::random_base()));
  mainp(string, sw_type<127>(tau, fp_type<127>::random_base()));

  return 0;
}
