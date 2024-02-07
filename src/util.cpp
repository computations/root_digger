#include <fstream>
#include <sstream>
#include <thread>
#include <unordered_set>
#if (!defined(__aarch64__))
#include <cpuid.h>
#endif
#ifdef MPI_BUILD
#include <mpi.h>
#endif
#include "util.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <sys/stat.h>

#ifdef _OPENMP
bool sysutil_dir_exists(const std::string &dname) {
  struct stat info;

  if (stat(dname.c_str(), &info) != 0) return false;
  else if (info.st_mode & S_IFDIR)
    return true;
  else
    return false;
}

std::string build_path(size_t cpu_number) {
  return "/sys/devices/system/cpu/cpu" + std::to_string(cpu_number)
         + "/topology/";
}

void get_cpuid(int32_t out[4], int32_t x) {
#ifdef _WIN32
  __cpuid(out, x);
#else
  __cpuid_count(x, 0, out[0], out[1], out[2], out[3]);
#endif
}

size_t read_id_from_file(const std::string &filename) {
  std::ifstream f(filename);
  if (f.good()) {
    size_t id;
    f >> id;
    return id;
  } else
    throw std::runtime_error("couldn't open sys files");
}

size_t get_numa_node_id(const std::string &cpu_path) {
  // this is ugly, but should be reliable -> please blame Linux kernel
  // developers & Intel!
  std::string node_path = cpu_path + "../node";
  for (size_t i = 0; i < 1000; ++i) {
    if (sysutil_dir_exists(node_path + std::to_string(i))) return i;
  }

  // fallback solution: return socket_id which is often identical to numa id
  return read_id_from_file(cpu_path + "physical_package_id");
}

size_t get_core_id(const std::string &cpu_path) {
  return read_id_from_file(cpu_path + "core_id");
}

size_t get_physical_core_count(size_t n_cpu) {
#if defined(__linux__)
  std::unordered_set<size_t> cores;
  for (size_t i = 0; i < n_cpu; ++i) {
    std::string cpu_path     = build_path(i);
    size_t      core_id      = get_core_id(cpu_path);
    size_t      node_id      = get_numa_node_id(cpu_path);
    size_t      uniq_core_id = (node_id << 16) + core_id;
    cores.insert(uniq_core_id);
  }
  return cores.size();
#else
  (void)(n_cpu);
  throw std::runtime_error("This function only supports linux");
#endif
}

bool ht_enabled() {
  int32_t info[4];

  get_cpuid(info, 1);

  return (bool)(info[3] & (0x1 << 28));
}

size_t sysutil_get_cpu_cores() {
  auto lcores = std::thread::hardware_concurrency();
  try {
    return get_physical_core_count(lcores);
  } catch (const std::runtime_error &) {
    auto threads_per_core = ht_enabled() ? 2u : 1u;

    return lcores / threads_per_core;
  }
}

#else
size_t sysutil_get_cpu_cores() { return 1; }
#endif

std::string combine_argv_argc(int argv, char **argc) {

  std::ostringstream oss;
  for (int i = 0; i < argv; ++i) {
    oss << argc[i];
    if (i != argv - 1) { oss << " "; }
  }
  return oss.str();
}
