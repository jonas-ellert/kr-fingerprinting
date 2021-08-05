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

template <MersennePrime p>
void mainp(int argc, char *argv[]) {
  static std::random_device seed;
  static std::mt19937_64 g(10);
  static std::uniform_int_distribution<uint8_t> d(1, 255);

  std::vector<uint8_t> string;

  if (argc > 1) {
    std::ifstream t(argv[1]);

    t.seekg(0, std::ios::end);
    string = std::vector<uint8_t>();
    string.reserve(t.tellg());
    t.seekg(0, std::ios::beg);

    string.assign((std::istreambuf_iterator<char>(t)),
                  std::istreambuf_iterator<char>());
    std::cout << "String loaded: " << string.size() << std::endl;
  } else {
    uint64_t const n = 1ULL * 512 * 1024 * 1024;
    string.resize(n);
    for (size_t i = 0; i < n; ++i) {
      string[i] = d(g);
    }
    std::cout << "String generated." << std::endl;
  }
  uint64_t const n = string.size();
  uint64_t const tau = std::stoi(argv[0]);

  [[maybe_unused]] auto set = [&](uint64_t const j, std::string const &s) {
    for (size_t i = 0; i < s.size(); i++) {
      string[i + j] = s[i];
    }
  };

  auto const b = kr_fingerprinting::kr_fingerprinter<p>::random_base();
  typename kr_fingerprinter<p>::sliding_window w(tau, b);
  std::cout << "b=" << b << std::endl;

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

  std::string name = std::string("FP") + std::to_string(p::shift);
  auto fp = measure(name, [&]() {
    typename p::uintX_t fp = 0;
    for (size_t i = 0; i < tau; i++) {
      fp = w.roll_right(fp, string[i]);
    }
    for (size_t i = 0; i < n - tau; i++) {
      fp = w.roll_right(fp, string[i], string[i + tau]);
    }
    return fp;
  });

  typename p::uintX_t fptest = 0;
  for (size_t i = n - tau; i < n; i++) {
    fptest = w.roll_right(fptest, string[i]);
  }
  std::cout << name << " correct=" << (fptest == fp) << std::endl;
}

int main(int argc, char *argv[]) {
  if (argc > 2 && std::string(argv[1]) == std::string("61")) {
    mainp<MERSENNE61>(argc - 2, &argv[2]);
  } else if (argc > 1 && std::string(argv[1]) == std::string("127")) {
    mainp<MERSENNE127>(argc - 2, &argv[2]);
  } else {
    mainp<MERSENNE61>(argc - 1, &argv[1]);
    mainp<MERSENNE89>(argc - 1, &argv[1]);
    mainp<MERSENNE107>(argc - 1, &argv[1]);
    mainp<MERSENNE127>(argc - 1, &argv[1]);
  }

  return 0;
}
