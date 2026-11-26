// dwarf_sycl.hpp -- shared SYCL host support for the sycl/ tree.
// Counterpart of opencl/include/common_args.c: same leading CLI ("-p P -d D -t T --"), device selection with a
// hard guard, ahead-of-region kernel build, and in-order-queue helpers that mirror blocking OpenCL enqueues.
#ifndef DWARF_SYCL_HPP
#define DWARF_SYCL_HPP
#include <sycl/sycl.hpp>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "dwarf_launch.h"
#include "dwarf_timer.h"

namespace dsycl {

struct Options { int platform = 0, device = -1, type = 0; };
inline Options g_opts;

// Same contract as ocd_parse(): nothing happens without "--"; otherwise argv advances so that argv[0] == "--".
inline void parse_args(int *argc, char ***argv) {
  int n = *argc, dd = -1;
  for (int i = 0; i < n; ++i)
    if (!std::strcmp((*argv)[i], "--")) { dd = i; break; }
  if (dd < 0) return;
  for (int i = 1; i < dd; ++i) {
    const char *a = (*argv)[i];
    if (!std::strcmp(a, "-o")) { std::fprintf(stderr, "-o (optimized OpenCL kernel variants) is not available in the SYCL tree\n"); std::exit(1); }
    if ((!std::strcmp(a, "-p") || !std::strcmp(a, "-d") || !std::strcmp(a, "-t")) && i + 1 < dd) {
      int v = std::atoi((*argv)[++i]);
      if (a[1] == 'p') g_opts.platform = v; else if (a[1] == 'd') g_opts.device = v; else g_opts.type = v;
    } else { std::fprintf(stderr, "unknown device option %s\nUsage: -p <platform> -d <device> -t <0 cpu|1 gpu> --\n", a); std::exit(1); }
  }
  *argv += dd; *argc -= dd;
}

inline void async_handler(sycl::exception_list el) {
  for (auto &e : el) {
    try { std::rethrow_exception(e); }
    catch (const std::exception &x) { std::fprintf(stderr, "SYCL asynchronous error: %s\n", x.what()); std::exit(1); }
  }
}

// Requested device type: DWARF_SYCL_DEVICE=cpu|gpu, otherwise -t (1 = gpu, anything else = cpu, as in OpenCL).
// -d indexes the devices of that type; -p is accepted for CLI compatibility and ignored.
inline sycl::queue make_queue() {
  const char *env = std::getenv("DWARF_SYCL_DEVICE");
  std::string want = (env && *env) ? env : (g_opts.type == 1 ? "gpu" : "cpu");
  if (want != "cpu" && want != "gpu") { std::fprintf(stderr, "DWARF_SYCL_DEVICE must be cpu or gpu\n"); std::exit(97); }
  std::vector<sycl::device> devs;
  try { devs = sycl::device::get_devices(want == "gpu" ? sycl::info::device_type::gpu : sycl::info::device_type::cpu); }
  catch (const sycl::exception &) {}
  int idx = g_opts.device < 0 ? 0 : g_opts.device;
  if (idx >= (int)devs.size()) { std::fprintf(stderr, "requested %s device %d not available (%zu found)\n", want.c_str(), idx, devs.size()); std::exit(97); }
  sycl::device dev = devs[idx];
  if ((want == "gpu") != dev.is_gpu() || (want == "cpu") != dev.is_cpu()) { std::fprintf(stderr, "selected device is not the requested %s\n", want.c_str()); std::exit(97); }
  std::string plat = dev.get_platform().get_info<sycl::info::platform::name>();
  std::printf("Platform Chosen : %s\n", plat.c_str());
  std::printf("Number of available devices: %zu\n", devs.size());
  std::printf("Device Chosen : %s\n", dev.get_info<sycl::info::device::name>().c_str());
  std::printf("SYCL requested=%s driver=%s\n", want.c_str(), dev.get_info<sycl::info::device::driver_version>().c_str());
  std::fflush(stdout);
  setenv("DWARF_DEVICE", want.c_str(), 1);   // label = the device actually in use
  setenv("DWARF_TOOLCHAIN", "dpcpp", 0);
  return sycl::queue(dev, async_handler, sycl::property::queue::in_order{});
}

// Counterpart of clBuildProgram: builds every kernel of the executable before the timed region.
// The bundle is kept alive in a static so that the runtime's program cache stays populated.
inline void prebuild(sycl::queue &q) {
  static std::vector<sycl::kernel_bundle<sycl::bundle_state::executable>> keep;
  keep.push_back(sycl::get_kernel_bundle<sycl::bundle_state::executable>(q.get_context(), {q.get_device()}));
}

template <class T> inline T *alloc(sycl::queue &q, size_t n) {
  T *p = sycl::malloc_device<T>(n ? n : 1, q);
  if (!p) { std::fprintf(stderr, "device allocation of %zu bytes failed\n", n * sizeof(T)); std::exit(1); }
  return p;
}
inline void release(sycl::queue &q, void *p) { sycl::free(p, q); }

// Blocking host->device copy (clEnqueueWriteBuffer with CL_TRUE).
inline void write(sycl::queue &q, void *dev, const void *host, size_t bytes) { q.memcpy(dev, host, bytes).wait_and_throw(); }

// When DWARF_DUMP_DIR is set, every device->host copy is appended to $DWARF_DUMP_DIR/readback.bin as
// {uint64 bytes, data}; sycl/tools/ocl_dump.c produces the same stream from the OpenCL binaries.
inline void dump(const void *host, size_t bytes) {
  static const char *dir = std::getenv("DWARF_DUMP_DIR");
  if (!dir || !*dir) return;
  static std::FILE *f = std::fopen((std::string(dir) + "/readback.bin").c_str(), "wb");
  if (!f) { std::perror("DWARF_DUMP_DIR"); std::exit(1); }
  unsigned long long n = bytes;
  std::fwrite(&n, sizeof n, 1, f); std::fwrite(host, 1, bytes, f); std::fflush(f);
}
// Blocking device->host copy (clEnqueueReadBuffer with CL_TRUE).
inline void read(sycl::queue &q, void *host, const void *dev, size_t bytes) {
  q.memcpy(host, dev, bytes).wait_and_throw();
  dump(host, bytes);
}

} // namespace dsycl
#endif
