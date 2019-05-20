#ifndef __RD_MODEL_HPP_
#define __RD_MODEL_HPP_

extern "C" {
#include <libpll/pll.h>
}
#include "tree.hpp"
#include <string>
#include <vector>

typedef std::vector<double> model_params_t;

std::string read_file_contents(std::ifstream &infile);

double parse_param(std::string::const_iterator begin,
                   std::string::const_iterator end);

model_params_t parse_model_params(const std::string &model_string);

model_params_t parse_model_file(const std::string &model_filename);

class model_t {
  model_t(const model_params_t &rate_parameters, pll_utree_t *tree,
          pll_msa_t *msa, unsigned int states);
  ~model_t();

private:
  pll_utree_t *_tree;
  pll_partition_t *_partition;
};

#endif
