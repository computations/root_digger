#ifndef __RD_MODEL_HPP_
#define __RD_MODEL_HPP_

extern "C" {
#include <libpll/pll.h>
}
#include "msa.hpp"
#include "tree.hpp"
#include <string>
#include <utility>
#include <vector>

typedef std::vector<double> model_params_t;

std::string read_file_contents(std::ifstream &infile);

double parse_param(std::string::const_iterator begin,
                   std::string::const_iterator end);

model_params_t parse_model_params(const std::string &model_string);

model_params_t parse_model_file(const std::string &model_filename);

class model_t {
public:
  model_t(const model_params_t &, rooted_tree_t, const msa_t &);
  ~model_t();
  double compute_lh(const root_location_t &root_location);
  double compute_lh_root(const root_location_t &root);
  double compute_dlh(const root_location_t &root_location);
  root_location_t optimize_alpha(const root_location_t &root);
  root_location_t optimize_root_location();
  const rooted_tree_t &rooted_tree(const root_location_t &root);

private:
  std::pair<root_location_t, double>
  bisect(const root_location_t &beg, double d_beg, const root_location_t &end,
         double d_end, double atol, size_t depth);

  rooted_tree_t _tree;
  pll_partition_t *_partition;
};

#endif
