#include "model.hpp"
#include <fstream>
#include <unordered_map>
extern "C" {
#include <libpll/pll_msa.h>
}

std::string read_file_contents(std::ifstream &infile) {
  std::string str;
  infile.seekg(0, std::ios::end);
  str.resize(infile.tellg());
  infile.seekg(0, std::ios::beg);
  infile.read(&str[0], str.size());
  return str;
}

double parse_param(std::string::const_iterator begin,
                   std::string::const_iterator end) {
  try {
    std::string tmp{begin, end};
    return std::stod(tmp);
  } catch (const std::invalid_argument &) {
    throw std::invalid_argument("Failed to parse model parameter string");
  } catch (const std::out_of_range &) {
    throw std::out_of_range(
        "Model parameter string contains a value which is too large");
  }
}

model_params_t parse_model_params(const std::string &model_string) {
  model_params_t model_params;
  auto substr_begin = model_string.begin();

  for (auto substr_end = model_string.begin();
       substr_begin != model_string.end(); ++substr_end) {

    if (*substr_end == ',') {
      model_params.push_back(parse_param(substr_begin, substr_end));
      substr_begin = substr_end + 1;
    } else if (*substr_end == '\0') {
      model_params.push_back(parse_param(substr_begin, substr_end));
      break;
    }
  }

  return model_params;
}

model_params_t parse_model_file(const std::string &model_filename) {
  std::ifstream model_file(model_filename);

  return parse_model_params(read_file_contents(model_file));
}

model_t::model_t(const model_params_t &rate_parameters, rooted_tree_t tree,
                 const msa_t &msa) {

  _tree = std::move(tree);

  /*
   * Only one submodel will be used for the time being. If there is desire for
   * more, we can add support for more models..
   */
  constexpr unsigned int submodels = 1;

  /*
   * For now, nonrev only supports the "generic" cpu case. I might be able to
   * boost it up to the other simd architectures for CLVS, but that is untested
   * at this moment.
   */
  unsigned int attributes = 0;
  attributes |= PLL_ATTRIB_ARCH_CPU;
  attributes |= PLL_ATTRIB_NONREV;

  _partition = pll_partition_create(_tree.tip_count(), _tree.branch_count(),
                                    msa.states(), msa.length(), submodels,
                                    _tree.branch_count(), submodels,
                                    _tree.branch_count(), attributes);
  pll_set_subst_params(_partition, 0, rate_parameters.data());

  /* set the frequencies
   *
   * For now we are going to just use emperical frequencies, but in the future,
   * we can do something more clever.
   */
  // double *emp_freqs = pllmod_msa_empirical_frequencies(_partition);
  std::vector<double> freqs{.25, .25, .25, .25};
  pll_set_frequencies(_partition, 0, freqs.data());
  // free(emp_freqs);

  double rate_cats[submodels] = {0};
  pll_compute_gamma_cats(1, 4, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates(_partition, rate_cats);

  /* update the invariant sites */
  pll_update_invariant_sites(_partition);

  /* no need to set rates and rate weights, only 1 category */

  /* make a label map */
  auto label_map = _tree.label_map();

  /* use the label map to assign tip states in the partition */
  for (int i = 0; i < msa.count(); ++i) {
    pll_set_tip_states(_partition, label_map.at(msa.label(i)), msa.map(),
                       msa.sequence(i));
  }

  /* set pattern weights */
  pll_set_pattern_weights(_partition, msa.weights());
}

model_t::~model_t() { pll_partition_destroy(_partition); }

double model_t::compute_lh(const root_location_t &root_location) {
  std::vector<pll_operation_t> ops;
  std::vector<unsigned int> pmatrix_indices;
  std::vector<double> branch_lengths;

  GENERATE_AND_UNPACK_OPS(_tree, root_location, ops, pmatrix_indices,
                          branch_lengths);

  unsigned int params[] = {0};

  int result =
      pll_update_prob_matrices(_partition, params, pmatrix_indices.data(),
                               branch_lengths.data(), pmatrix_indices.size());

  if (result == PLL_FAILURE) {
    throw std::runtime_error(pll_errmsg);
  }

  pll_update_partials(_partition, ops.data(), ops.size());

  double lh = pll_compute_root_loglikelihood(_partition, _tree.root_clv_index(),
                                             _tree.root_scaler_index(), params,
                                             nullptr);
  return lh;
}
