#include "model.hpp"
#include <cmath>
#include <fstream>
#include <string>
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

/*
 * Use a tangent method to compute the derivative
 */
double model_t::compute_dlh(const root_location_t &root) {
  constexpr double EPSILON = 1e-4;

  root_location_t root_prime{root};
  root_prime.brlen_ratio += EPSILON;
  if (root_prime.brlen_ratio >= 1.0) {
    throw std::runtime_error(
        "Alpha out of bounds while computing derivative: " +
        std::to_string(root.brlen_ratio));
  }

  double fx = compute_lh(root);
  if (!std::isfinite(fx)) {
    throw std::runtime_error("fx is not finite when computing derivative: " +
                             std::to_string(root.edge->length));
  }
  double fxh = compute_lh(root_prime);
  if (!std::isfinite(fxh)) {
    throw std::runtime_error("fxh is not finite when computing derivative: " +
                             std::to_string(root_prime.edge->length));
  }

  return (fxh - fx) / EPSILON;
}

std::pair<root_location_t, double> model_t::bisect(const root_location_t &beg,
                                                   double d_beg,
                                                   const root_location_t &end,
                                                   double d_end, double atol,
                                                   size_t depth = 0) {
  root_location_t midpoint{beg};
  midpoint.brlen_ratio = (beg.brlen_ratio + end.brlen_ratio) / 2;

  double d_midpoint = compute_dlh(midpoint);

  if (depth > 64) {
    return {midpoint, compute_lh(midpoint)};
  }

  /* case 1: d_midpoint is within tolerance of 0.0
   * We found a root, and should return
   */

  if (fabs(d_midpoint) < atol) {
    return {midpoint, compute_lh(midpoint)};
  }

  /* case 2: midpoint is opposite sign of both beg and end
   * There are 2 roots, so we recurse on both sides, and return the one with the
   * best LH
   */
  if ((d_beg > 0.0 && d_end > 0.0 && d_midpoint < 0.0) ||
      (d_beg < 0.0 && d_end < 0.0 && d_midpoint > 0.0)) {
    auto r1 = bisect(beg, d_beg, midpoint, d_midpoint, atol, depth + 1);
    auto r2 = bisect(midpoint, d_midpoint, end, d_end, atol, depth + 1);
    if (r1.second < r2.second) {
      return r2;
    }
    return r1;
  }
  /* case 3: end and midpoint share a sign, while beg has the opposite sign
   * In this case, there is a root between midpoint and beg
   */
  if ((d_beg < 0.0 && d_midpoint > 0.0 && d_end > 0.0) ||
      (d_beg > 0.0 && d_midpoint < 0.0 && d_end < 0.0)) {
    return bisect(beg, d_beg, midpoint, d_midpoint, atol, depth + 1);
  }

  /* case 4: beg and midpoint share a sign, while end has the opposite sign
   * In this case, there is a root between midpoint and end
   */
  if ((d_beg < 0.0 && d_midpoint < 0.0 && d_end > 0.0) ||
      (d_beg > 0.0 && d_midpoint > 0.0 && d_end < 0.0)) {
    return bisect(midpoint, d_midpoint, end, d_end, atol, depth + 1);
  }

  /* case 5: something went wrong */
  throw std::runtime_error(
      "Bisection failed to converge with interval : " + std::to_string(d_beg) +
      ", " + std::to_string(d_end) +
      ", midpoint: " + std::to_string(d_midpoint));
}

/* Find the optimum for the ratio via the bisection method */
root_location_t model_t::optimize_alpha(const root_location_t &root) {
  constexpr double ATOL = 1e-6;
  root_location_t beg{root};
  beg.brlen_ratio = 1e-4;

  root_location_t end{root};
  end.brlen_ratio = 1.0 - 1e-4 * 2;

  double d_beg = compute_dlh(beg);
  if (fabs(d_beg) < ATOL) {
    return beg;
  }

  double d_end = compute_dlh(end);
  if (fabs(d_end) < ATOL) {
    return end;
  }

  if (!std::isfinite(d_beg) || !std::isfinite(d_end)) {
    throw std::runtime_error(
        "Initial derivatives failed when optimizing alpha: " +
        std::to_string(root.edge->length));
  }

  if ((d_beg < 0.0 && d_end > 0.0) || (d_beg > 0.0 && d_end < 0.0)) {
    return bisect(beg, d_beg, end, d_end, ATOL).first;
  }

  /* if we have gotten here, then we must have an "even" function so, grid
   * search to find a midpoint which has an opposite sign value and use that to
   * do bisection with
   */

  bool beg_end_pos = d_beg > 0.0 && d_end > 0.0;

  for (size_t midpoints = 2; midpoints <= 16; midpoints *= 2) {
    for (size_t midpoint = 1; midpoint <= midpoints; ++midpoint) {
      if (midpoint % 2 == 0)
        continue;

      double alpha = 1.0 / (double)midpoints * midpoint;
      root_location_t midpoint_root{beg};
      midpoint_root.brlen_ratio = alpha;
      double d_midpoint = compute_dlh(midpoint_root);
      if (fabs(d_midpoint) < ATOL) {
        return midpoint_root;
      }
      if ((beg_end_pos && d_midpoint < 0.0) ||
          (!beg_end_pos && d_midpoint > 0.0)) {
        auto r1 = bisect(beg, d_beg, midpoint_root, d_midpoint, ATOL);
        auto r2 = bisect(midpoint_root, d_midpoint, end, d_end, ATOL);
        if (r1.second < r2.second) {
          return r2.first;
        }
        return r1.first;
      }
    }
  }

  throw std::runtime_error(
      "Initial derivatives failed when optimizing alpha, ran out of cases");
}
