extern "C" {
#include <libpll/pll.h>
}
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
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

void print_usage() {
  std::cout << "Usage: rd [options]\n"
            << "Version: " << GIT_REV_STRING << "\n"
            << "Application Options:\n"
            << "    --msa [FILE]\n"
            << "           File containing the alignment.\n"
            << "    --tree [NUMBER]\n"
            << "           File containing the tree, with branch lengths.\n"
            << "    --model [RATIO]\n"
            << "           File containing the substitution model.\n"
            << "           Ideally, this should be a non-reversible model.\n"
            << "    --freqs [FILE]\n"
            << "           Optional file containing the frequencies. If not "
               "given, then\n"
            << "           empirical frequencies are used instead.\n"
            << "    --verbose\n"
            << "       Enable debug output. Warning, extremely noisy\n"
            << std::endl;
}

int main(int argv, char **argc) {
  static struct option long_opts[] = {
      {"msa", required_argument, 0, 0},
      {"tree", required_argument, 0, 0},
      {"model", required_argument, 0, 0},
      {"freqs", required_argument, 0, 0},
      {"verbose", no_argument, 0, 0},
      {"version", no_argument, 0, 0},
      {0, 0, 0, 0},
  };

  try {
    int c;
    int index = 0;
    std::string msa_filename;
    std::string tree_filename;
    std::string model_filename;
    std::string freqs_filename;
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
      case 4: // verbose
        __VERBOSE__ = true;
        break;
      case 5: //version
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
    }
    if (model_filename.empty()) {
      throw std::invalid_argument("A model file is required");
    }

    model_params_t params = parse_model_file(model_filename);
    model_params_t freqs;
    if (!freqs_filename.empty()) {
      freqs = parse_model_file(freqs_filename);
    }

    unsigned int states =
        static_cast<unsigned int>(std::ceil(std::sqrt(params.size())));
    const pll_state_t *map = nullptr;

    if (states == 4)
      map = pll_map_nt;
    if (states == 20)
      map = pll_map_aa;

    if (map == nullptr) {
      throw std::invalid_argument(
          "Root digger only supports protein and nucleotide data");
    }

    msa_t msa{msa_filename, map, states};
    rooted_tree_t tree{tree_filename};
    if (freqs.size() == 0) {
      model_t model{params, tree, msa};
      auto rl = model.optimize_root_location();
      std::cout << model.rooted_tree(rl).newick() << std::endl;
    } else {
      model_t model{params, tree, msa, freqs};
      auto rl = model.optimize_root_location();
      std::cout << model.rooted_tree(rl).newick() << std::endl;
    }

  } catch (const std::exception &e) {
    std::cout << "There was an error during processing:\n"
              << e.what() << std::endl;
    return 1;
  }
  return 0;
}
