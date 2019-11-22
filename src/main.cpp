extern "C" {
#include <libpll/pll.h>
}
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "debug.h"
int __VERBOSE__ = EMIT_LEVEL_WARNING;
#include "model.hpp"
#include "msa.hpp"
#include "tree.hpp"

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)
#define GIT_REV_STRING STRINGIFY(GIT_REV)
#define GIT_COMMIT_STRING STRINGIFY(GIT_COMMIT)
#define BUILD_DATE_STRING STRINGIFY(BUILD_DATE)

void print_version() {
  /*
  std::cout << "Version: " << GIT_REV_STRING << "\n"
            << "Build Commit: " << GIT_COMMIT_STRING << "\n"
            << "Build Date: " << BUILD_DATE_STRING << std::endl;
  */
  debug_print(EMIT_LEVEL_IMPORTANT, "Version: %s", GIT_REV_STRING);
  debug_print(EMIT_LEVEL_IMPORTANT, "Build Commit: %s", GIT_COMMIT_STRING);
  debug_print(EMIT_LEVEL_IMPORTANT, "Build Date: %s", BUILD_DATE_STRING);
}

void print_run_header(
    const std::chrono::time_point<std::chrono::system_clock> &start_time,
    uint64_t seed) {
  time_t st = std::chrono::system_clock::to_time_t(start_time);
  char time_string[256];
  std::strftime(time_string, sizeof(time_string), "%F %T", std::localtime(&st));
  // std::cout << "Runing Root Digger \n";
  debug_string(EMIT_LEVEL_IMPORTANT, "Running Root Digger");
  print_version();
  /*
  std::cout << "Started: " << std::put_time(std::localtime(&st), "%F %T")
            << "\n"
            << "Seed: " << seed << std::endl;
            */
  debug_print(EMIT_LEVEL_IMPORTANT, "Started: %s", time_string);
  debug_print(EMIT_LEVEL_IMPORTANT, "Seed: %lu", seed);
}

void print_usage() {
  std::cout
      << "Usage: rd [options]\n"
      << "Version: " << GIT_REV_STRING << "\n"
      << "Application Options:\n"
      << "    --msa [FILE]\n"
      << "           File containing the alignment.\n"
      << "    --tree [FILE]\n"
      << "           File containing the tree, with branch lengths.\n"
      << "    --model [FILE]\n"
      << "           Optional file containing the substitution model.\n"
      << "           Ideally, this should be a non-reversible model.\n"
      << "           If this is not given, then the parameters are inferred\n"
      << "    --freqs [FILE]\n"
      << "           Optional file containing the frequencies. If not\n"
      << "           given, then empirical frequencies are used instead.\n"
      << "    --partition [FILE]\n"
      << "           Optional file containing the partition specification.\n"
      << "           Format is the same as RAxML-NG partition file.\n"
      << "    --seed [NUMBER]\n"
      << "           Random seed to use. Optional\n"
      << "    --min-roots [NUMBER]\n"
      << "           Minimum number of roots to start from. Optional,\n"
      << "           Default is 1.\n"
      << "    --root-ratio [NUMBER]\n"
      << "           Proportion of potential starting roots to attempt\n"
      << "           Default is 0.0\n"
      << "    --atol [NUMBER]\n"
      << "           Root optmization stopping tolerance. Increase this to \n"
      << "           improve results.Default is 1e-4\n"
      << "    --silent\n"
      << "           Suppress output except for the final tree\n"
      << "    --force\n"
      << "           Disable Saftey checks.\n"
      << "    --verbose\n"
      << "           Increase the verbosity level. Can be repeated to\n"
      << "           level further.\n"
      << std::endl;
}

int main(int argv, char **argc) {
  auto start_time = std::chrono::system_clock::now();
  static struct option long_opts[] = {
      {"msa", required_argument, 0, 0},        /* 0 */
      {"tree", required_argument, 0, 0},       /* 1 */
      {"model", required_argument, 0, 0},      /* 2 */
      {"freqs", required_argument, 0, 0},      /* 3 */
      {"seed", required_argument, 0, 0},       /* 4 */
      {"states", required_argument, 0, 0},     /* 5 */
      {"verbose", no_argument, 0, 0},          /* 6 */
      {"silent", no_argument, 0, 0},           /* 7 */
      {"min-roots", required_argument, 0, 0},  /* 8 */
      {"root-ratio", required_argument, 0, 0}, /* 9 */
      {"atol", required_argument, 0, 0},       /* 10 */
      {"bfgstol", required_argument, 0, 0},    /* 11 */
      {"factor", required_argument, 0, 0},     /* 12 */
      {"partition", required_argument, 0, 0},  /* 13 */
      {"exhaustive", no_argument, 0, 0},       /* 14 */
      {"force", no_argument, 0, 0},            /* 15 */
      {"version", no_argument, 0, 0},          /* 16 */
      {"debug", no_argument, 0, 0},            /* 17 */
      {0, 0, 0, 0},
  };

  if (argv == 1) {
    print_usage();
    return 0;
  }
  try {
    int c;
    int index = 0;
    std::string msa_filename;
    std::string tree_filename;
    std::string model_filename;
    std::string freqs_filename;
    std::string partition_filename;
    uint64_t seed = std::random_device()();
    size_t min_roots = 1;
    double root_ratio = 0.05;
    double abs_tolerance = 1e-7;
    double factor = 1e4;
    double bfgs_tol = 1e-7;
    unsigned int states = 0;
    bool silent = false;
    bool sanity_checks = true;
    bool exhaustive = false;
    while ((c = getopt_long_only(argv, argc, "", long_opts, &index)) == 0) {
      debug_print(EMIT_LEVEL_DEBUG, "parsing option index: %d", index);
      switch (index) {
      case 0: // msa
        msa_filename = optarg;
        break;
      case 1: // tree
        tree_filename = optarg;
        break;
      case 2: // model
        model_filename = optarg;
        break;
      case 3: // freqs
        freqs_filename = optarg;
        break;
      case 4: // seed
        seed = atol(optarg);
        break;
      case 5: // states
        states = atol(optarg);
        break;
      case 6: // verbose
        __VERBOSE__ += 1;
        break;
      case 7: // silent
        __VERBOSE__ = 0;
        silent = true;
        break;
      case 8: // min-roots
        min_roots = atol(optarg);
        break;
      case 9: // root-ratio
        root_ratio = atof(optarg);
        break;
      case 10: // atol
        abs_tolerance = atof(optarg);
        break;
      case 11: // bfgs_tol
        bfgs_tol = atof(optarg);
        break;
      case 12: // factor
        factor = atof(optarg);
        break;
      case 13: // partition
        partition_filename = optarg;
        break;
      case 14: // force
        exhaustive = true;
        break;
      case 15: // force
        sanity_checks = false;
        break;
      case 16: // version
        print_version();
        return 0;
      case 17: // debug
        __VERBOSE__ = EMIT_LEVEL_DEBUG;
        break;
      default:
        throw std::invalid_argument("An argument was not recognized");
      }
    }

    if (msa_filename.empty()) {
      print_usage();
      return 1;
    }
    if (tree_filename.empty()) {
      print_usage();
      return 1;
    }

    model_params_t model_params;

    if (!model_filename.empty()) {
      model_params = parse_model_file(model_filename);
      states =
          static_cast<unsigned int>(std::ceil(std::sqrt(model_params.size())));
    } else if (states == 0) {
      std::cout << "Please give either a model file, or a number of states to "
                   "the model"
                << std::endl;
      print_usage();
      return 1;
    }

    if (!silent)
      print_run_header(start_time, seed);

    debug_print(EMIT_LEVEL_INFO, "abs_tolerance: %f", abs_tolerance);
    const pll_state_t *map = nullptr;

    if (states == 4)
      map = pll_map_nt;
    if (states == 20)
      map = pll_map_aa;

    if (map == nullptr) {
      throw std::invalid_argument(
          "Root digger only supports protein and nucleotide data");
    }

    std::vector<msa_t> msa;
    if (partition_filename.empty()) {
      msa.emplace_back(msa_filename, map, states);
    } else {
      msa_t unparted_msa{msa_filename, map, states};
      msa = unparted_msa.partition(parse_partition_file(partition_filename));
    }

    rooted_tree_t tree{tree_filename};

    if (sanity_checks) {
      if (!tree.sanity_check()) {
        std::cout << "WARNING: Refusing to run. An extremely long branch is "
                     "present in the provided tree. Roots are extremely likely "
                     "to be erronously placed on this branch. Use '--force' to "
                     "override this warning"
                  << std::endl;
        return 1;
      }
    }

    model_t model{tree, msa, seed};
    try {
      model.initialize_partitions(msa);
    } catch (const invalid_empirical_frequencies_exception &) {
      model.initialize_partitions_uniform_freqs(msa);
    }
    root_location_t final_rl;
    double final_lh;
    if (!exhaustive) {
      auto tmp = model.optimize_all(min_roots, root_ratio, abs_tolerance,
                                    bfgs_tol, factor);
      final_rl = tmp.first;
      final_lh = tmp.second;
    } else {
      auto tmp = model.exhaustive_search(abs_tolerance, bfgs_tol, factor);
      final_rl = tmp.first;
      final_lh = tmp.second;
    }
    std::string final_tree_string{model.rooted_tree(final_rl).newick()};
    if (!silent) {
      debug_print(EMIT_LEVEL_IMPORTANT, "Final LogLH: %.5f", final_lh);
    }
    std::cout << final_tree_string << std::endl;
    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    if (!silent)
      std::cout << "Inference took: " << duration.count() << "s" << std::endl;

  } catch (const std::exception &e) {
    std::cout << "There was an error during processing:\n"
              << e.what() << std::endl;
    return 1;
  }

  return 0;
}
