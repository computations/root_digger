extern "C" {
#include <libpll/pll.h>
}
#ifndef _WIN32
#include <cpuid.h>
#endif
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#ifdef MPI_BUILD
#include <mpi.h>
#endif
#include <omp.h>
#include <random>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <vector>

#include "debug.h"
int __VERBOSE__ = EMIT_LEVEL_PROGRESS;
int __MPI_RANK__ = 0;
int __MPI_NUM_TASKS__ = 1;
#include "model.hpp"
#include "msa.hpp"
#include "tree.hpp"

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)
#define GIT_REV_STRING STRINGIFY(GIT_REV)
#define GIT_COMMIT_STRING STRINGIFY(GIT_COMMIT)
#define BUILD_DATE_STRING STRINGIFY(BUILD_DATE)

static void print_version() {
  debug_print(EMIT_LEVEL_IMPORTANT, "Version: %s", GIT_REV_STRING);
  debug_print(EMIT_LEVEL_IMPORTANT, "Build Commit: %s", GIT_COMMIT_STRING);
  debug_print(EMIT_LEVEL_IMPORTANT, "Build Date: %s", BUILD_DATE_STRING);
}

static std::string combine_argv_argc(int argv, char **argc) {

  std::ostringstream oss;
  for (int i = 0; i < argv; ++i) {
    oss << argc[i];
    if (i != argv - 1) {
      oss << " ";
    }
  }
  return oss.str();
}

static void print_run_header(
    const std::chrono::time_point<std::chrono::system_clock> &start_time,
    uint64_t seed, size_t threads, int argv, char **argc) {
  time_t st = std::chrono::system_clock::to_time_t(start_time);
  char time_string[256];
  std::strftime(time_string, sizeof(time_string), "%F %T", std::localtime(&st));
  debug_string(EMIT_LEVEL_IMPORTANT, "Running Root Digger");
  print_version();
  debug_print(EMIT_LEVEL_IMPORTANT, "Started: %s", time_string);
  debug_print(EMIT_LEVEL_IMPORTANT, "Seed: %lu", seed);
  debug_print(EMIT_LEVEL_IMPORTANT, "Number of threads per proc: %lu", threads);
#ifdef MPI_VERSION
  debug_print(EMIT_LEVEL_IMPORTANT, "Number of procs %d", __MPI_NUM_TASKS__);
#endif
  debug_print(EMIT_LEVEL_IMPORTANT, "Command: %s",
              combine_argv_argc(argv, argc).c_str());
}

class initialized_flag_t {
public:
  enum value_t {
    uninitalized,
    initialized_true,
    initialized_false,
  };
  initialized_flag_t() : value(value_t::uninitalized){};
  constexpr initialized_flag_t(const value_t &v) : value(v) {}
  constexpr bool operator==(const initialized_flag_t &rhs) const {
    return rhs.value == value;
  }
  constexpr bool operator!=(const initialized_flag_t &rhs) const {
    return rhs.value != value;
  }
  constexpr bool initalized() const { return value != value_t::uninitalized; }
  initialized_flag_t &operator=(const initialized_flag_t &rhs) {
    value = rhs.value;
    return *this;
  }
  bool convert_with_default(bool default_value) {
    if (value == value_t::uninitalized)
      return default_value;
    return value == value_t::initialized_true;
  }

private:
  value_t value;
};

/* The following functions are taken from RAxML-NG. I did some of the initial
 * work, but the final implementation is from Alexey Kozlov
 */

#ifdef _OPENMP
bool sysutil_dir_exists(const std::string &dname) {
  struct stat info;

  if (stat(dname.c_str(), &info) != 0)
    return false;
  else if (info.st_mode & S_IFDIR)
    return true;
  else
    return false;
}

static std::string build_path(size_t cpu_number) {
  return "/sys/devices/system/cpu/cpu" + std::to_string(cpu_number) +
         "/topology/";
}

static void get_cpuid(int32_t out[4], int32_t x) {
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
    if (sysutil_dir_exists(node_path + std::to_string(i)))
      return i;
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
    std::string cpu_path = build_path(i);
    size_t core_id = get_core_id(cpu_path);
    size_t node_id = get_numa_node_id(cpu_path);
    size_t uniq_core_id = (node_id << 16) + core_id;
    cores.insert(uniq_core_id);
  }
  return cores.size();
#else
  (void)(n_cpu);
  throw std::runtime_error("This function only supports linux");
#endif
}

static bool ht_enabled() {
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

struct cli_options_t {
  std::string msa_filename;
  std::string tree_filename;
  std::string output_tree_filename;
  std::string model_filename;
  std::string freqs_filename;
  std::string partition_filename;
  std::string data_type;
  std::string model_string;
  std::vector<size_t> rate_cats{1};
  std::vector<rate_category::rate_category_e> rate_category_types;
  uint64_t seed = std::random_device()();
  size_t min_roots = 1;
  size_t threads = 0;
  double root_ratio = 0.01;
  double abs_tolerance = 1e-7;
  double factor = 1e4;
  double br_tolerance = 1e-12;
  double bfgs_tol = 1e-7;
  const unsigned int states = 4;
  bool silent = false;
  bool exhaustive = false;
  bool echo = false;
  bool invariant_sites = false;
  initialized_flag_t early_stop;
};

static void print_usage() {
  if (__MPI_RANK__ != 0) {
    return;
  }
  std::cout
      << "Usage: rd [options]\n"
      << "Version: " << GIT_REV_STRING << "\n"
      << "Application Options:\n"
      << "    --msa [FILE]\n"
      << "           File containing the alignment.\n"
      << "    --tree [FILE]\n"
      << "           File containing the tree, with branch lengths.\n"
      << "    --partition [FILE]\n"
      << "           Optional file containing the partition specification.\n"
      << "           Format is the same as RAxML-NG partition file.\n"
      << "    --exhaustive\n"
      << "           Enable exhaustive mode. This will attempt to root a tree\n"
      << "           at every branch, and then report the results using LWR.\n"
      << "    --early-stop\n"
      << "           Enable early stopping. This will cause cause the search\n"
      << "           to terminate when the root placement is sufficently\n"
      << "           close for 2 consecutive iterations. How close they need\n"
      << "           to be is controled by brtol. Is enabled by default for\n"
      << "           search mode and disabled by default for exhaustive mode.\n"
      << "    --no-early-stop\n"
      << "           Force disable early stop.\n"
      << "    --seed [NUMBER]\n"
      << "           Random seed to use. Optional\n"
      << "    --rate-cats [NUMBER]\n"
      << "           Number of rate categories to use for the model. Default\n"
      << "           is 1.\n"
      << "    --invariant-sites\n"
      << "           Enable invariant sites. Default is off.\n"
      << "    --min-roots [NUMBER]\n"
      << "           Minimum number of roots to start from. Optional,\n"
      << "           Default is 1.\n"
      << "    --root-ratio [NUMBER]\n"
      << "           Proportion of potential starting roots to attempt\n"
      << "           Default is 0.01\n"
      << "    --atol [NUMBER]\n"
      << "           Root optmization stopping tolerance. Increase this to \n"
      << "           improve results.Default is 1e-4\n"
      << "    --brtol [NUMBER]\n"
      << "           When early stop mode is enabled, this controls the\n"
      << "           distance required to trigger. Default is 1e-12\n"
      << "    --bfgstol [NUMBER]\n"
      << "           Tolerance for the BFGS steps. Default is 1e-7\n"
      << "    --factor [NUMBER]\n"
      << "           Factor for the BFGS steps. Default is 1e4\n"
      << "    --threads [NUMBER]\n"
      << "           Number of threads to use\n"
      << "    --silent\n"
      << "           Suppress output except for the final tree\n"
      << "    --verbose\n"
      << "           Increase the verbosity level. Can be repeated to\n"
      << "           level further.\n"
      << std::endl;
}

int wrapped_main(int argv, char **argc) {
#ifdef MPI_VERSION
  MPI_Comm_size(MPI_COMM_WORLD, &__MPI_NUM_TASKS__);
  MPI_Comm_rank(MPI_COMM_WORLD, &__MPI_RANK__);
#endif

  auto start_time = std::chrono::system_clock::now();
  static struct option long_opts[] = {
      {"msa", required_argument, 0, 0},             /* 0 */
      {"tree", required_argument, 0, 0},            /* 1 */
      {"model", required_argument, 0, 0},           /* 2 */
      {"seed", required_argument, 0, 0},            /* 3 */
      {"verbose", no_argument, 0, 0},               /* 4 */
      {"silent", no_argument, 0, 0},                /* 5 */
      {"min-roots", required_argument, 0, 0},       /* 6 */
      {"root-ratio", required_argument, 0, 0},      /* 7 */
      {"atol", required_argument, 0, 0},            /* 8 */
      {"brtol", required_argument, 0, 0},           /* 9 */
      {"bfgstol", required_argument, 0, 0},         /* 10 */
      {"factor", required_argument, 0, 0},          /* 11 */
      {"partition", required_argument, 0, 0},       /* 12 */
      {"treefile", required_argument, 0, 0},        /* 13 */
      {"exhaustive", no_argument, 0, 0},            /* 14 */
      {"early-stop", no_argument, 0, 0},            /* 15 */
      {"no-early-stop", no_argument, 0, 0},         /* 16 */
      {"rate-cats", required_argument, 0, 0},       /* 17 */
      {"rate-cats-type", required_argument, 0, 0},  /* 18 */
      {"invariant-sites", required_argument, 0, 0}, /* 19 */
      {"threads", required_argument, 0, 0},         /* 20 */
      {"version", no_argument, 0, 0},               /* 21 */
      {"debug", no_argument, 0, 0},                 /* 22 */
      {"mpi-debug", no_argument, 0, 0},             /* 22 */
      {"echo", no_argument, 0, 0},                  /* 23 */
      {0, 0, 0, 0},
  };

  if (argv == 1) {
    print_usage();
    return 0;
  }
  try {
    int c;
    int index = 0;
    cli_options_t cli_options;
    while ((c = getopt_long_only(argv, argc, "", long_opts, &index)) == 0) {
      debug_print(EMIT_LEVEL_DEBUG, "parsing option index: %d", index);
      switch (index) {
      case 0: // msa
        cli_options.msa_filename = optarg;
        break;
      case 1: // tree
        cli_options.tree_filename = optarg;
        break;
      case 2: // model
        cli_options.model_string = optarg;
        break;
      case 3: // seed
        cli_options.seed = static_cast<uint64_t>(atol(optarg));
        break;
      case 4: // verbose
        __VERBOSE__ += 1;
        break;
      case 5: // silent
        __VERBOSE__ = 0;
        cli_options.silent = true;
        break;
      case 6: // min-roots
        cli_options.min_roots = static_cast<size_t>(atol(optarg));
        break;
      case 7: // root-ratio
        cli_options.root_ratio = atof(optarg);
        break;
      case 8: // atol
        cli_options.abs_tolerance = atof(optarg);
        break;
      case 9: // brtol
        cli_options.br_tolerance = atof(optarg);
        break;
      case 10: // bfgs_tol
        cli_options.bfgs_tol = atof(optarg);
        break;
      case 11: // factor
        cli_options.factor = atof(optarg);
        break;
      case 12: // partition
        cli_options.partition_filename = optarg;
        break;
      case 13: // treefile
        cli_options.output_tree_filename = optarg;
        break;
      case 14: // exhaustive
        cli_options.exhaustive = true;
        break;
      case 15: // early-stop
        cli_options.early_stop = initialized_flag_t::initialized_true;
        break;
      case 16: // no-early-stop
        cli_options.early_stop = initialized_flag_t::initialized_false;
        break;
      case 17: // rate-cats
        cli_options.rate_cats = {(size_t)atol(optarg)};
        break;
      case 18: // rate-cats-type
        if (strcmp(optarg, "mean") == 0) {
          cli_options.rate_category_types.push_back(rate_category::MEAN);
        } else if (strcmp(optarg, "median") == 0) {
          cli_options.rate_category_types.push_back(rate_category::MEDIAN);
        } else if (strcmp(optarg, "free") == 0) {
          cli_options.rate_category_types.push_back(rate_category::FREE);
        }
        break;
      case 19: // invariant-sites
        cli_options.invariant_sites = true;
        break;
      case 20: // threads
        cli_options.threads = {(size_t)atol(optarg)};
        break;
      case 21: // version
        print_version();
        return 0;
      case 22: // debug
        __VERBOSE__ = EMIT_LEVEL_DEBUG;
        break;
      case 23: // mpi-debug
        __VERBOSE__ = EMIT_LEVEL_MPI_DEBUG;
        break;
      case 24: // echo
        cli_options.echo = true;
        break;
      default:
        throw std::invalid_argument("An argument was not recognized");
      }
    }
#ifdef MPI_VERSION
    debug_string(EMIT_LEVEL_DEBUG, "Broadcasting seed");
    MPI_Bcast(&cli_options.seed, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    debug_print(EMIT_LEVEL_MPI_DEBUG, "seed: %lu", cli_options.seed);
#endif

    if (cli_options.msa_filename.empty()) {
      std::cout << "No MSA was given, please supply an MSA" << std::endl;
      print_usage();
      return 1;
    }
    if (cli_options.tree_filename.empty()) {
      std::cout << "No tree was given, please supply an tree" << std::endl;
      print_usage();
      return 1;
    }

    if (cli_options.threads == 0) {
#ifdef MPI_VERSION
      cli_options.threads = 1;
#else
      cli_options.threads = sysutil_get_cpu_cores();
#endif
    }

#ifdef _OPENMP
    omp_set_num_threads(cli_options.threads);
#endif

    if (!cli_options.silent)
      print_run_header(start_time, cli_options.seed, cli_options.threads, argv,
                       argc);

    debug_print(EMIT_LEVEL_INFO, "abs_tolerance: %.08f",
                cli_options.abs_tolerance);
    if (cli_options.exhaustive &&
        cli_options.early_stop.convert_with_default(!cli_options.exhaustive)) {
    }

    const pll_state_t *map = nullptr;
    if (cli_options.states == 4)
      map = pll_map_nt;
    else if (cli_options.states == 20)
      map = pll_map_aa;
    else
      throw std::runtime_error("Data Type is not supported");

    if (map == nullptr) {
      throw std::invalid_argument(
          "Root digger only supports protein and nucleotide data");
    }

    if (!cli_options.model_string.empty()) {
      auto mi = parse_model_info(cli_options.model_string);
      cli_options.rate_cats.clear();
      cli_options.rate_category_types.clear();
      cli_options.rate_cats.push_back(mi.ratehet_opts.rate_cats);
      cli_options.rate_category_types.push_back(
          mi.ratehet_opts.rate_category_type);
      if (mi.ratehet_opts.alpha_init) {
        debug_string(EMIT_LEVEL_WARNING,
                     "Ignoring alpha in model string as it currently "
                     "is not suported");
      }
      auto subst_str{mi.subst_str};

      for (auto &ch : subst_str) {
        ch = std::tolower(ch);
      }
      if (subst_str != "unrest") {
        debug_print(EMIT_LEVEL_WARNING,
                    "Ignoring subst matrix %s for model from command line"
                    ". Currently only UNREST is supported",
                    mi.subst_str.c_str());
      }
    }

    std::vector<msa_t> msa;
    msa_partitions_t part_infos;
    if (cli_options.partition_filename.empty()) {
      msa.emplace_back(cli_options.msa_filename, map, cli_options.states);
    } else {
      msa_t unparted_msa{cli_options.msa_filename, map, cli_options.states,
                         false};
      part_infos = parse_partition_file(cli_options.partition_filename);
      msa = unparted_msa.partition(part_infos);
    }

    if (part_infos.size() > 0) {
      if (cli_options.rate_cats.size() > 0) {
        debug_string(EMIT_LEVEL_WARNING,
                     "Using rate categories from the partition file "
                     "over the option passed on the command line");
      }
      cli_options.rate_cats.clear();
      for (auto &p : part_infos) {

        size_t rate_cats = p.model.ratehet_opts.rate_cats;
        if (rate_cats == 0) {
          rate_cats = 1;
        }
        cli_options.rate_cats.push_back(rate_cats);
        if (p.model.ratehet_opts.alpha_init) {
          debug_print(EMIT_LEVEL_WARNING,
                      "Ignoring alpha in partition %s as it currently "
                      "is not suported",
                      p.partition_name.c_str());
        }
        auto subst_str{p.model.subst_str};

        for (auto &ch : subst_str) {
          ch = std::tolower(ch);
        }

        if (subst_str != "unrest") {
          debug_print(EMIT_LEVEL_WARNING,
                      "Ignoring subst matrix %s for partition "
                      "%s. Currently only UNREST is supported",
                      p.model.subst_str.c_str(), p.partition_name.c_str());
        }
      }
    }

    if (cli_options.rate_category_types.size() == 0) {
      cli_options.rate_category_types.push_back(rate_category::MEAN);
    }

    if (cli_options.rate_cats.size() == 1 && cli_options.rate_cats[0] == 0) {
      throw std::runtime_error("Rate categories cannot be zero");
    }

    for (auto &m : msa) {
      m.valid_data();
    }

    rooted_tree_t tree{cli_options.tree_filename};

    if (cli_options.min_roots > tree.root_count()) {
      throw std::runtime_error(
          "Min roots is larger than the number of roots on the tree");
    }

    if (cli_options.root_ratio < 0) {
      throw std::runtime_error("Root ratio is negative");
    }

#ifdef MPI_VERSION
    if (__MPI_NUM_TASKS__ == 1) {
      debug_string(EMIT_LEVEL_WARNING,
                   "Running MPI version with only 1 process, "
                   "is this really what you meant?");
    }
#endif

    model_t model{
        tree,
        msa,
        cli_options.rate_cats,
        cli_options.rate_category_types,
        cli_options.invariant_sites,
        cli_options.seed,
        cli_options.early_stop.convert_with_default(!cli_options.exhaustive)};
    try {
      model.initialize_partitions(msa);
    } catch (const invalid_empirical_frequencies_exception &) {
      model.initialize_partitions_uniform_freqs(msa);
    }

    if (cli_options.echo) {
      std::cout << tree.newick() << std::endl;
    }

    model.compute_lh(tree.root_location(0));
    root_location_t final_rl;
    double final_lh = -std::numeric_limits<double>::infinity();
    std::string final_tree_string;
    if (!cli_options.exhaustive) {
      model.assign_indicies_by_rank_search(
          cli_options.min_roots, cli_options.root_ratio,
          static_cast<size_t>(__MPI_RANK__),
          static_cast<size_t>(__MPI_NUM_TASKS__));
      auto tmp =
          model.optimize_all(cli_options.min_roots, cli_options.root_ratio,
                             cli_options.abs_tolerance, cli_options.bfgs_tol,
                             cli_options.br_tolerance, cli_options.factor);
      if (__MPI_RANK__ == 0) {
        final_rl = tmp.first;
        final_lh = tmp.second;
        final_tree_string = model.rooted_tree(final_rl).newick();
      }
    } else {

      model.assign_indicies_by_rank_exhaustive(
          static_cast<size_t>(__MPI_RANK__),
          static_cast<size_t>(__MPI_NUM_TASKS__));

      auto tmp = model.exhaustive_search(
          cli_options.abs_tolerance, cli_options.bfgs_tol,
          cli_options.br_tolerance, cli_options.factor);
      if (__MPI_RANK__ == 0) {
        final_rl = tmp.first;
        final_lh = tmp.second;
        final_tree_string = model.virtual_rooted_tree(final_rl).newick();
      }
    }
    if (!cli_options.silent) {
      debug_print(EMIT_LEVEL_IMPORTANT, "Final LogLH: %.5f", final_lh);
    }
    if (__MPI_RANK__ == 0) {
      std::cout << final_tree_string << std::endl;
    }
    auto end_time = std::chrono::system_clock::now();

    std::chrono::duration<double> duration = end_time - start_time;

    if (!cli_options.silent && __MPI_RANK__ == 0)
      std::cout << "Inference took: " << duration.count() << "s" << std::endl;

    if (!cli_options.output_tree_filename.empty() && __MPI_RANK__ == 0) {
      std::ofstream outfile{cli_options.output_tree_filename};
      outfile << final_tree_string;
    }

  } catch (const std::exception &e) {
    if (__MPI_RANK__ == 0) {
      std::cout << "There was an error during processing:\n"
                << e.what() << std::endl;
    }
    return 1;
  }

  return 0;
}

int main(int argv, char **argc) {
#ifdef MPI_VERSION
  MPI_Init(&argv, &argc);
#endif
  wrapped_main(argv, argc);
#ifdef MPI_VERSION
  MPI_Finalize();
#endif
}
