extern "C" {
#include <libpll/pll.h>
}
#include <chrono>
#include <cmath>
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
#include "model.hpp"
#include "msa.hpp"
#include "tree.hpp"

bool __VERBOSE__ = false;

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)
#define GIT_REV_STRING STRINGIFY(GIT_REV)
#define GIT_COMMIT_STRING STRINGIFY(GIT_COMMIT)
#define BUILD_DATE_STRING STRINGIFY(BUILD_DATE)

void print_version() {
  std::cout << "Version: " << GIT_REV_STRING << "\n"
            << "Build Commit: " << GIT_COMMIT_STRING << "\n"
            << "Build Date: " << BUILD_DATE_STRING << std::endl;
}

void print_run_header(
    const std::chrono::time_point<std::chrono::system_clock> &start_time,
    uint64_t seed) {
  time_t st = std::chrono::system_clock::to_time_t(start_time);
  std::cout << "Runing Root Digger \n";
  print_version();
  std::cout << "Started: " << std::put_time(std::localtime(&st), "%F %T")
            << "\n"
            << "Seed: " << seed << std::endl;
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
      << "    --silent\n"
      << "           Suppress output except for the final tree\n"
      << "    --fast\n"
      << "           Use a fast annealing schedule. Gives results quickly,\n"
      << "           but might be less optimal.\n"
      << "    --slow\n"
      << "           Use a slow annealing schedule. Gives results slowly, \n"
      << "           but results are more likely to be optimial.\n"
      << "    --tfinal [NUMBER]\n"
      << "           Threshold to end the simulated annealing at.\n"
      << "           Default is 1e-6\n"
      << "    --force\n"
      << "           Disable Saftey checks.\n"
      << "    --verbose\n"
      << "           Enable debug output. Warning, extremely noisy\n"
      << std::endl;
}

int main(int argv, char **argc) {
  auto start_time = std::chrono::system_clock::now();
  static struct option long_opts[] = {
      {"msa", required_argument, 0, 0},
      {"tree", required_argument, 0, 0},
      {"model", required_argument, 0, 0},
      {"freqs", required_argument, 0, 0},
      {"seed", required_argument, 0, 0},
      {"states", required_argument, 0, 0},
      {"verbose", no_argument, 0, 0},
      {"silent", no_argument, 0, 0},
      {"fast", no_argument, 0, 0},
      {"slow", no_argument, 0, 0},
      {"force", no_argument, 0, 0},
      {"tfinal", no_argument, 0, 0},
      {"partition", required_argument, 0, 0},
      {"root-freq", required_argument, 0, 0},
      {"version", no_argument, 0, 0},
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
    unsigned int states = 0;
    bool silent = false;
    double temp_param = 0.9;
    double root_opt_freq = 2.0;
    bool sanity_checks = true;
    while ((c = getopt_long_only(argv, argc, "", long_opts, &index)) == 0) {
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
        __VERBOSE__ = true;
        break;
      case 7: // silent
        silent = true;
        break;
      case 8: // fast
        temp_param = 0.5;
        break;
      case 9: // slow
        temp_param = 0.98;
        break;
      case 10: // force
        sanity_checks = false;
        break;
      case 11: // final_temp
        // final_temp = atof(optarg);
        break;
      case 12: // partition
        partition_filename = optarg;
        break;
      case 13: // root-freq
        root_opt_freq = atof(optarg);
        break;
      case 14: // version
        print_version();
        return 0;
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

    std::string final_tree_string;
    double final_lh = -INFINITY;

    for (size_t i = 0; i < 1; ++i) {
      model_t model{tree, msa, seed + i};
      try {
        model.initialize_partitions(msa);
      } catch (const invalid_empirical_frequencies_exception &) {
        model.initialize_partitions_uniform_freqs(msa);
      }
      model.set_temp_ratio(temp_param);
      model.set_root_opt_frequency(root_opt_freq);
      auto rl = model.optimize_all();
      double lh = model.compute_lh(rl);
      if (lh > final_lh) {
        final_tree_string = model.rooted_tree(rl).newick();
        final_lh = lh;
      }
    }
    std::cout << std::fixed << std::setprecision(2);
    std::cout << final_lh << std::endl;
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
