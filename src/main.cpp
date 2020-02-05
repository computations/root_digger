extern "C" {
#include <libpll/pll.h>
}
#include <algorithm>
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
int __VERBOSE__ = EMIT_LEVEL_PROGRESS;
#include "model.hpp"
#include "msa.hpp"
#include "tree.hpp"

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)
#define GIT_REV_STRING STRINGIFY(GIT_REV)
#define GIT_COMMIT_STRING STRINGIFY(GIT_COMMIT)
#define BUILD_DATE_STRING STRINGIFY(BUILD_DATE)

void print_version() {
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
  debug_string(EMIT_LEVEL_IMPORTANT, "Running Root Digger");
  print_version();
  debug_print(EMIT_LEVEL_IMPORTANT, "Started: %s", time_string);
  debug_print(EMIT_LEVEL_IMPORTANT, "Seed: %lu", seed);
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
  constexpr bool operator==(const initialized_flag_t &rhs) {
    return rhs.value == value;
  }
  constexpr bool operator!=(const initialized_flag_t &rhs) {
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

struct cli_options_t {
  std::string msa_filename;
  std::string tree_filename;
  std::string output_tree_filename;
  std::string model_filename;
  std::string freqs_filename;
  std::string partition_filename;
  std::string data_type;
  std::vector<size_t> rate_cats{1};
  uint64_t seed = std::random_device()();
  size_t min_roots = 1;
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

void print_usage() {
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
      << "    --silent\n"
      << "           Suppress output except for the final tree\n"
      << "    --verbose\n"
      << "           Increase the verbosity level. Can be repeated to\n"
      << "           level further.\n"
      << std::endl;
}

int main(int argv, char **argc) {
  auto start_time = std::chrono::system_clock::now();
  static struct option long_opts[] = {
      {"msa", required_argument, 0, 0},             /* 0 */
      {"tree", required_argument, 0, 0},            /* 1 */
      {"model", required_argument, 0, 0},           /* 2 */
      {"seed", required_argument, 0, 0},            /* 4 */
      {"verbose", no_argument, 0, 0},               /* 5 */
      {"silent", no_argument, 0, 0},                /* 6 */
      {"min-roots", required_argument, 0, 0},       /* 7 */
      {"root-ratio", required_argument, 0, 0},      /* 8 */
      {"atol", required_argument, 0, 0},            /* 9 */
      {"brtol", required_argument, 0, 0},           /* 10 */
      {"bfgstol", required_argument, 0, 0},         /* 11 */
      {"factor", required_argument, 0, 0},          /* 12 */
      {"partition", required_argument, 0, 0},       /* 13 */
      {"treefile", required_argument, 0, 0},        /* 14 */
      {"exhaustive", no_argument, 0, 0},            /* 15 */
      {"early-stop", no_argument, 0, 0},            /* 16 */
      {"no-early-stop", no_argument, 0, 0},         /* 17 */
      {"rate-cats", required_argument, 0, 0},       /* 18 */
      {"invariant-sites", required_argument, 0, 0}, /* 19 */
      {"version", no_argument, 0, 0},               /* 20 */
      {"debug", no_argument, 0, 0},                 /* 21 */
      {"echo", no_argument, 0, 0},                  /* 22 */
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
        cli_options.model_filename = optarg;
        break;
      case 3: // seed
        cli_options.seed = atol(optarg);
        break;
      case 4: // verbose
        __VERBOSE__ += 1;
        break;
      case 5: // silent
        __VERBOSE__ = 0;
        cli_options.silent = true;
        break;
      case 6: // min-roots
        cli_options.min_roots = atol(optarg);
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
      case 18:
        cli_options.invariant_sites = true;
        break;
      case 19: // version
        print_version();
        return 0;
      case 20: // debug
        __VERBOSE__ = EMIT_LEVEL_DEBUG;
        break;
      case 21: // echo
        cli_options.echo = true;
        break;
      default:
        throw std::invalid_argument("An argument was not recognized");
      }
    }

    if (cli_options.msa_filename.empty()) {
      print_usage();
      return 1;
    }
    if (cli_options.tree_filename.empty()) {
      print_usage();
      return 1;
    }

    model_params_t model_params;

    if (!cli_options.model_filename.empty()) {
      model_params = parse_model_file(cli_options.model_filename);
    } else if (cli_options.states == 0) {
      std::cout << "Please give either a model file, or a number of states to "
                   "the model"
                << std::endl;
      print_usage();
      return 1;
    }

    if (!cli_options.silent)
      print_run_header(start_time, cli_options.seed);

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

    std::vector<msa_t> msa;
    if (cli_options.partition_filename.empty()) {
      msa.emplace_back(cli_options.msa_filename, map, cli_options.states);
    } else {
      msa_t unparted_msa{cli_options.msa_filename, map, cli_options.states};
      msa = unparted_msa.partition(
          parse_partition_file(cli_options.partition_filename));
    }

    rooted_tree_t tree{cli_options.tree_filename};

    model_t model{
        tree,
        msa,
        cli_options.rate_cats,
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

    root_location_t final_rl;
    double final_lh;
    std::string final_tree_string;
    if (!cli_options.exhaustive) {
      auto tmp =
          model.optimize_all(cli_options.min_roots, cli_options.root_ratio,
                             cli_options.abs_tolerance, cli_options.bfgs_tol,
                             cli_options.br_tolerance, cli_options.factor);
      final_rl = tmp.first;
      final_lh = tmp.second;
      final_tree_string = std::move(model.rooted_tree(final_rl).newick());
    } else {
      auto tmp = model.exhaustive_search(
          cli_options.abs_tolerance, cli_options.bfgs_tol,
          cli_options.br_tolerance, cli_options.factor);
      final_rl = tmp.first;
      final_lh = tmp.second;
      final_tree_string =
          std::move(model.virtual_rooted_tree(final_rl).newick());
    }
    if (!cli_options.silent) {
      debug_print(EMIT_LEVEL_IMPORTANT, "Final LogLH: %.5f", final_lh);
    }
    std::cout << final_tree_string << std::endl;
    auto end_time = std::chrono::system_clock::now();

    std::chrono::duration<double> duration = end_time - start_time;

    if (!cli_options.silent)
      std::cout << "Inference took: " << duration.count() << "s" << std::endl;

    if (!cli_options.output_tree_filename.empty()) {
      std::ofstream outfile{cli_options.output_tree_filename};
      outfile << final_tree_string;
    }

  } catch (const std::exception &e) {
    std::cout << "There was an error during processing:\n"
              << e.what() << std::endl;
    return 1;
  }

  return 0;
}
