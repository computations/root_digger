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

struct dlh_t {
  double lh;
  double dlh;
};

std::string read_file_contents(std::ifstream &infile);

double parse_param(std::string::const_iterator begin,
                   std::string::const_iterator end);

model_params_t parse_model_params(const std::string &model_string);

model_params_t parse_model_file(const std::string &model_filename);

model_params_t random_params(size_t size, uint64_t seed);

class model_t {
public:
  model_t(model_params_t, rooted_tree_t, const msa_t &, const model_params_t &,
          uint64_t);
  model_t(model_params_t p, rooted_tree_t t, const msa_t &msa, uint64_t seed)
      : model_t(p, t, msa, model_params_t(), seed){};
  model_t(rooted_tree_t t, const msa_t &msa, model_params_t f, uint64_t seed)
      : model_t(std::move(random_params(
                    msa.states() * msa.states() - msa.states(), seed)),
                t, msa, f, seed){};
  model_t(rooted_tree_t t, const msa_t &msa, uint64_t seed)
      : model_t(std::move(random_params(
                    msa.states() * msa.states() - msa.states(), seed)),
                t, msa, model_params_t(), seed) {}
  ~model_t();
  double compute_lh(const root_location_t &root_location);
  double compute_lh_root(const root_location_t &root);
  dlh_t compute_dlh(const root_location_t &root_location);
  root_location_t optimize_alpha(const root_location_t &root);
  std::pair<root_location_t, double> optimize_root_location();
  root_location_t optimize_all(double final_temp);
  const rooted_tree_t &rooted_tree(const root_location_t &root);


  void set_temp_ratio(double);
  std::string subst_string() const;

private:
  std::pair<root_location_t, double>
  bisect(const root_location_t &beg, dlh_t d_beg, const root_location_t &end,
         dlh_t d_end, double atol, size_t depth);

  model_params_t _subst_params;
  rooted_tree_t _tree;
  pll_partition_t *_partition;
  uint64_t _seed;
  double _temp_ratio;
};

#endif
