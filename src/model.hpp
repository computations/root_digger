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

struct invalid_emperical_frequencies_exception : public std::runtime_error {
  invalid_emperical_frequencies_exception(const char *m)
      : std::runtime_error(m){};
};

std::string read_file_contents(std::ifstream &infile);

double parse_param(std::string::const_iterator begin,
                   std::string::const_iterator end);

model_params_t parse_model_params(const std::string &model_string);

model_params_t parse_model_file(const std::string &model_filename);

model_params_t random_params(size_t size, uint64_t seed);

class model_t {
public:
  model_t(rooted_tree_t t, const std::vector<msa_t> &msa, uint64_t seed);
  ~model_t();
  double compute_lh(const root_location_t &root_location);
  double compute_lh_root(const root_location_t &root);
  dlh_t compute_dlh(const root_location_t &root_location);
  root_location_t optimize_alpha(const root_location_t &root);
  std::pair<root_location_t, double> optimize_root_location();
  root_location_t optimize_all(double final_temp);
  const rooted_tree_t &rooted_tree(const root_location_t &root);

  void initialize_partitions(const std::vector<msa_t> &);
  void initialize_partitions_uniform_freqs(const std::vector<msa_t> &);
  void set_temp_ratio(double);
  std::string subst_string() const;

private:
  std::pair<root_location_t, double>
  bisect(const root_location_t &beg, dlh_t d_beg, const root_location_t &end,
         dlh_t d_end, double atol, size_t depth);
  std::pair<root_location_t, double> brents(root_location_t beg, dlh_t d_beg,
                                            root_location_t end, dlh_t d_end,
                                            double atol);
  void set_subst_rates(size_t, const model_params_t &);
  void set_subst_rates_random(size_t, const msa_t &);
  void set_gamma_rates(size_t);
  void update_invariant_sites(size_t);
  void set_tip_states(size_t, const msa_t &);
  void set_empirical_freqs(size_t);
  void set_freqs(size_t, const model_params_t &);
  void anneal_rates(const std::vector<model_params_t> &,
                    const std::vector<model_params_t> &, 
                    const root_location_t&,double, double);

  std::vector<model_params_t> _subst_params;
  rooted_tree_t _tree;
  std::vector<pll_partition_t *> _partitions;
  std::vector<unsigned int> _partition_weights;
  uint64_t _seed;
  double _temp_ratio;
};

#endif
