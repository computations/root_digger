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

#include "model.hpp"
#include "msa.hpp"
#include "tree.hpp"

int main(int argv, char **argc) {
  static struct option long_opts[] = {
      {"msa", required_argument, 0, 0},
      {"tree", required_argument, 0, 0},
      {"model", required_argument, 0, 0},
      {0, 0, 0, 0},
  };

  int c;
  int index = 0;
  std::string msa_filename;
  std::string tree_filename;
  std::string model_filename;
  while ((c = getopt_long_only(argv, argc, "", long_opts, &index)) == 0) {
    switch (index) {
    case 0:
      msa_filename = optarg;
      break;
    case 1:
      tree_filename = optarg;
      break;
    case 2:
      model_filename = optarg;
      break;
    }
  }

  model_params_t params = parse_model_file(model_filename);

  unsigned int states = static_cast<unsigned int>(std::ceil(std::sqrt(params.size())));
  const pll_state_t *map = nullptr;

  if (states == 4)
    map = pll_map_nt;
  if (states == 20)
    map = pll_map_aa;

  if (map == nullptr) {
    throw std::invalid_argument(
        "Root digger only supports protien and nuculotide data");
  }

  msa_t msa{msa_filename, map, states};
  rooted_tree_t tree{tree_filename};
  model_t model{params, tree, msa};
  auto rl = model.optimize_root_location();
  std::cout << model.rooted_tree(rl).newick() << std::endl;

  return 0;
}
