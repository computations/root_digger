extern "C" {
#include <libpll/pll.h>
}
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

  char c;
  int index = 0;
  std::string msa_filename;
  std::string tree_filename;
  std::string model_filename;
  while ((c = getopt_long_only(argv, argc, "", long_opts, &index))) {
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

  return 0;
}
