#ifndef DWARFS_OMP_COMMON_DWARF_COMMON_HPP
#define DWARFS_OMP_COMMON_DWARF_COMMON_HPP

#include <algorithm>
#include <cerrno>
#include <climits>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifndef DWARF_THREAD_LIMIT
#define DWARF_THREAD_LIMIT 256
#endif

namespace dwarfs_omp {

struct Options {
  std::size_t size = 1024;
  std::size_t secondary_size = 0;
  int repeats = 3;
  int warmup = 1;
  int iterations = 1;
  bool verify = true;
  std::string input_file;
  std::string input_file2;
  std::string input_file3;
};

inline bool MatchArg(const char *arg, const char *name) {
  return std::strcmp(arg, name) == 0;
}

inline Options ParseOptions(int argc, char **argv) {
  Options opts;
  auto fail = [](const std::string &message) {
    std::cerr << "ERROR: " << message << "\nUse --help for options.\n";
    std::exit(2);
  };
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto value = [&]() -> const char * {
      if (i + 1 >= argc) fail("missing value for " + arg);
      return argv[++i];
    };
    auto number = [&](int minimum) -> int {
      const char *raw = value();
      char *end = nullptr;
      errno = 0;
      const long long n = std::strtoll(raw, &end, 10);
      if (errno || end == raw || *end || n < minimum || n > INT_MAX)
        fail("invalid value for " + arg + ": " + raw);
      return static_cast<int>(n);
    };
    if (arg == "--size" || arg == "-n") opts.size = number(1);
    else if (arg == "--secondary-size") opts.secondary_size = number(0);
    else if (arg == "--repeats") opts.repeats = number(1);
    else if (arg == "--warmup") opts.warmup = number(0);
    else if (arg == "--iterations") opts.iterations = number(1);
    else if (arg == "--no-verify") opts.verify = false;
    else if (arg == "--verify") opts.verify = true;
    else if (arg == "--input" || arg == "-i") opts.input_file = value();
    else if (arg == "--input2") opts.input_file2 = value();
    else if (arg == "--input3") opts.input_file3 = value();
    else if (arg == "--help" || arg == "-h") {
      std::cout << "Options: --size N --secondary-size N --iterations N "
                   "--warmup N --repeats N [--input FILE] [--input2 FILE] "
                   "[--input3 FILE] [--verify | --no-verify]\n";
      std::exit(0);
    } else fail("unknown option: " + arg);
  }
  if (opts.secondary_size == 0) opts.secondary_size = opts.size;
  return opts;
}

class Timer {
 public:
  void Start() { start_ = Clock::now(); }
  double StopMs() const {
    const auto end = Clock::now();
    return std::chrono::duration<double, std::milli>(end - start_).count();
  }

 private:
  using Clock = std::chrono::steady_clock;
  Clock::time_point start_;
};

// Opt-in phase-boundary logging for cross-process alignment (e.g. an
// external power-sampling process integrating only over real kernel
// windows, not whole-process wall time which is dominated by device-init
// at small sizes). Silent/zero-overhead unless DWARF_PHASE_LOG is set in
// the environment -- existing stdout (DWARF_RESULT) is never touched,
// phase lines go to stderr only, so nothing that parses stdout is affected.
// Uses system_clock (wall-clock/epoch, comparable across processes)
// deliberately, NOT Timer's steady_clock (process-local, no fixed origin,
// unusable for cross-process alignment) -- the two clocks serve different
// purposes and are not meant to be unified.
class PhaseLogger {
 public:
  static PhaseLogger &Instance() {
    static PhaseLogger instance;
    return instance;
  }

  void Mark(const char *phase, int rep, const char *event) {
    if (!enabled_) return;
    const auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    std::cerr << "DWARF_PHASE,phase=" << phase << ",rep=" << rep
             << ",event=" << event << ",wall_ns=" << now_ns << "\n";
  }

 private:
  PhaseLogger() : enabled_(std::getenv("DWARF_PHASE_LOG") != nullptr) {}
  bool enabled_;
};

// Run the requested warmups and measured repetitions.
// setup_fn: untimed preparation called before each rep (pass [](){} if none needed).
// run_fn:   the timed portion, including any inner opts.iterations loop.
template <typename SetupFn, typename RunFn>
inline std::vector<double> AdaptiveRun(const Options &opts,
                                        SetupFn setup_fn, RunFn run_fn) {
  auto &plog = PhaseLogger::Instance();
  // Warmup reps get their OWN setup/kernel marks (rep < 0, so they never
  // collide with the timed reps' rep=0..N-1) instead of one undifferentiated
  // "warmup" blob -- the first-ever kernel submission (real device-init/JIT
  // cost included) happens here, and a consumer that wants that cost
  // (device-performance timing per this project's phase-boundary convention:
  // first pragma region submitted to final cross-device synchronize, excluding
  // only host-side I/O) needs it separated from setup_fn()'s own time the
  // same way the timed reps below already are.
  for (int i = 0; i < opts.warmup; ++i) {
    const int rep = -(i + 1);
    plog.Mark("setup", rep, "begin");
    setup_fn();
    plog.Mark("setup", rep, "end");
    plog.Mark("kernel", rep, "begin");
    run_fn();
    plog.Mark("kernel", rep, "end");
  }
  plog.Mark("setup", 0, "begin");
  setup_fn();
  plog.Mark("setup", 0, "end");
  Timer probe;
  probe.Start();
  plog.Mark("kernel", 0, "begin");
  run_fn();
  plog.Mark("kernel", 0, "end");
  const double probe_ms = probe.StopMs();
  const int repeats = opts.repeats;
  std::vector<double> times;
  times.reserve(repeats);
  times.push_back(probe_ms);
  for (int r = 1; r < repeats; ++r) {
    plog.Mark("setup", r, "begin");
    setup_fn();
    plog.Mark("setup", r, "end");
    Timer t;
    t.Start();
    plog.Mark("kernel", r, "begin");
    run_fn();
    plog.Mark("kernel", r, "end");
    times.push_back(t.StopMs());
  }
  return times;
}

inline double Mean(const std::vector<double> &values) {
  double sum = 0.0;
  for (double v : values) {
    sum += v;
  }
  return values.empty() ? 0.0 : sum / static_cast<double>(values.size());
}

inline double Median(std::vector<double> values) {
  if (values.empty()) {
    return 0.0;
  }
  std::sort(values.begin(), values.end());
  const std::size_t mid = values.size() / 2;
  if (values.size() % 2 == 0) {
    return 0.5 * (values[mid - 1] + values[mid]);
  }
  return values[mid];
}

inline const char *TargetName() {
#if defined(DWARF_TARGET_NAME)
  return DWARF_TARGET_NAME;
#else
  return "unknown";
#endif
}

inline void PrintResult(const std::string &benchmark, const std::string &motif,
                        std::size_t size, std::size_t secondary_size,
                        int iterations, int repeats,
                        const std::vector<double> &times_ms, bool valid) {
  std::cout << "DWARF_RESULT"
            << ",benchmark=" << benchmark
            << ",motif=" << motif
            << ",target=" << TargetName()
            << ",size=" << size
            << ",secondary_size=" << secondary_size
            << ",iterations=" << iterations
            << ",repeats=" << repeats
            << ",mean_ms=" << Mean(times_ms)
            << ",median_ms=" << Median(times_ms)
            << ",min_ms=" << (times_ms.empty() ? 0.0 : *std::min_element(times_ms.begin(), times_ms.end()))
            << ",max_ms=" << (times_ms.empty() ? 0.0 : *std::max_element(times_ms.begin(), times_ms.end()))
            << ",valid=" << (valid ? 1 : 0) << "\n";
}

template <typename T>
inline bool NearlyEqual(T a, T b, double tolerance = 1.0e-3) {
  const double da = static_cast<double>(a);
  const double db = static_cast<double>(b);
  const double denom = std::max(1.0, std::max(std::abs(da), std::abs(db)));
  return std::abs(da - db) / denom <= tolerance;
}

}  // namespace dwarfs_omp

#endif  // DWARFS_OMP_COMMON_DWARF_COMMON_HPP
