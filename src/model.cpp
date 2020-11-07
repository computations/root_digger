#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "checkpoint.hpp"
#include "debug.h"
#include "libpll/pll.h"
#include "model.hpp"
#include "msa.hpp"
#include "pll.h"
#include "tree.hpp"
#include "util.hpp"
extern "C" {
#include <lbfgsb.h>
#include <libpll/pll_msa.h>
}

std::string read_file_contents(std::ifstream &infile) {
  std::string str;
  infile.seekg(0, std::ios::end);
  str.resize(static_cast<unsigned long>(infile.tellg()));
  infile.seekg(0, std::ios::beg);
  infile.read(&str[0], static_cast<long>(str.size()));
  return str;
}

static std::string to_string(const std::vector<double> &v) {
  std::stringstream out;
  out << std::setprecision(3);
  out << "[";
  for (auto &i : v) {
    out << std::setw(3);
    out << std::fixed << i << " ";
  }
  out.seekp(-1, out.cur);
  out << "]";
  return out.str();
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
  auto           substr_begin = model_string.begin();

  for (auto substr_end = model_string.begin();
       substr_begin != model_string.end();
       ++substr_end) {
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
  if (!model_file) {
    throw std::runtime_error{"Failed to open the model file"};
  }

  return parse_model_params(read_file_contents(model_file));
}

model_params_t random_params(size_t size, uint64_t seed) {
  model_params_t                   mp(size);
  std::minstd_rand                 engine(seed);
  std::uniform_real_distribution<> dist(1e-4, 1.0);
  for (auto &f : mp) { f = dist(engine); }
  return mp;
}

model_t::model_t(
    rooted_tree_t                                      tree,
    const std::vector<msa_t> &                         msas,
    const std::vector<size_t> &                        rate_cats,
    const std::vector<rate_category::rate_category_e> &rate_category_types,
    bool                                               invariant_sites,
    uint64_t                                           seed,
    bool                                               early_stop) :
    _rate_category_types{rate_category_types},
    _invariant_sites{invariant_sites},
    _seed{seed},
    _early_stop{early_stop} {
  if (_early_stop) {
    debug_string(EMIT_LEVEL_IMPORTANT, "INFO: Early stop is enabled");
  }

  _random_engine = std::minstd_rand(_seed);
  _tree          = std::move(tree);
  for (auto rc : rate_cats) {
    _rate_rates.emplace_back(rc, 1.0);
    _rate_weights.emplace_back(rc, 1.0 / rc);
    _param_indicies.emplace_back(rc, 0);
  }

  for (auto &msa : msas) {
    if (!msa.constiency_check(_tree.label_set())) {
      throw std::invalid_argument(
          "Taxa on the tree and in the MSA are inconsistient");
    }
  }

  unsigned int attributes = 0;
  if (PLL_STAT(avx2_present)) {
    attributes |= PLL_ATTRIB_ARCH_AVX2;
  } else if (PLL_STAT(avx_present)) {
    attributes |= PLL_ATTRIB_ARCH_AVX;
  } else if (PLL_STAT(sse42_present)) {
    attributes |= PLL_ATTRIB_ARCH_SSE;
  }

  attributes |= PLL_ATTRIB_NONREV;
  attributes |= PLL_ATTRIB_SITE_REPEATS;

  size_t total_weight = 0;
  for (size_t partition_index = 0; partition_index < msas.size();
       ++partition_index) {
    auto &msa = msas[partition_index];

    if (_rate_rates[partition_index].size()
        > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("The size of the rate weights for the partition "
                               "is too large to safely cast");
    }

    if (msa.length() > static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error(
          "The length of the MSA is too large to safely cast");
    }

    _partitions.push_back(pll_partition_create(
        _tree.tip_count(),
        _tree.branch_count(),
        msa.states(),
        msa.length(),
        _submodels,
        _tree.branch_count(),
        static_cast<unsigned int>(_rate_rates[partition_index].size()),
        _tree.branch_count(),
        attributes));
    _partition_weights.push_back(msa.total_weight());

    set_gamma_rates(partition_index);
    total_weight += msa.total_weight();
  }

  assign_indicies();
}

model_t::~model_t() {
  for (auto p : _partitions) {
    if (p) pll_partition_destroy(p);
  }
}

void model_t::set_subst_rates(size_t p_index, const model_params_t &mp) {
  pll_set_subst_params(_partitions[p_index], 0, mp.data());
}

void model_t::set_subst_rates_random(size_t p_index, const msa_t &msa) {
  set_subst_rates(p_index,
                  random_params(msa.states() * msa.states() - msa.states(),
                                _random_engine()));
}

void model_t::set_subst_rates_random(size_t p_index, size_t states) {
  set_subst_rates(p_index,
                  random_params(states * states - states, _random_engine()));
}

void model_t::set_gamma_weights(size_t p_index, model_params_t w) {
  double sum = 0.0;
  debug_print(EMIT_LEVEL_DEBUG, "input weights %s", to_string(w).c_str());
  for (auto &f : w) { sum += f; }
  for (auto &f : w) { f /= sum; }
  debug_print(EMIT_LEVEL_DEBUG, "setting weights to %s", to_string(w).c_str());
  pll_set_category_weights(_partitions[p_index], w.data());
}

void model_t::set_gamma_rates(size_t p_index) {
  pll_set_category_weights(_partitions[p_index], _rate_weights[p_index].data());
  switch (_rate_category_types[p_index]) {
  case rate_category::MEAN:
    set_gamma_rates_mean(p_index);
    break;
  case rate_category::MEDIAN:
    set_gamma_rates_median(p_index);
    break;
  default:
    set_gamma_rates_free(p_index);
    break;
  }
}

void model_t::set_gamma_rates(size_t p_index, const model_params_t &alpha) {
  switch (_rate_category_types[p_index]) {
  case rate_category::MEAN:
    set_gamma_rates_mean(p_index, alpha[0]);
    break;
  case rate_category::MEDIAN:
    set_gamma_rates_median(p_index, alpha[0]);
    break;
  default:
    set_gamma_rates_free(p_index, alpha);
    break;
  }
}

void model_t::set_gamma_rates_mean(size_t p_index) {
  pll_compute_gamma_cats(1.0,
                         static_cast<unsigned int>(_rate_rates[p_index].size()),
                         _rate_rates[p_index].data(),
                         PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates(_partitions[p_index], _rate_rates[p_index].data());
}

void model_t::set_gamma_rates_mean(size_t p_index, double alpha) {
  pll_compute_gamma_cats(alpha,
                         static_cast<unsigned int>(_rate_rates[p_index].size()),
                         _rate_rates[p_index].data(),
                         PLL_GAMMA_RATES_MEDIAN);
  pll_set_category_rates(_partitions[p_index], _rate_rates[p_index].data());
}

void model_t::set_gamma_rates_median(size_t p_index) {
  pll_compute_gamma_cats(1.0,
                         static_cast<unsigned int>(_rate_rates[p_index].size()),
                         _rate_rates[p_index].data(),
                         PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates(_partitions[p_index], _rate_rates[p_index].data());
}

void model_t::set_gamma_rates_median(size_t p_index, double alpha) {
  pll_compute_gamma_cats(alpha,
                         static_cast<unsigned int>(_rate_rates[p_index].size()),
                         _rate_rates[p_index].data(),
                         PLL_GAMMA_RATES_MEDIAN);
  pll_set_category_rates(_partitions[p_index], _rate_rates[p_index].data());
}

void model_t::set_gamma_rates_free(size_t p_index) {
  for (auto &r : _rate_rates[p_index]) { r = 1.0; }
  pll_set_category_rates(_partitions[p_index], _rate_rates[p_index].data());
}

void model_t::set_gamma_rates_free(size_t p_index, model_params_t free_rates) {
  double sum = 0.0;
  debug_print(
      EMIT_LEVEL_DEBUG, "input rates %s", to_string(free_rates).c_str());
  for (size_t i = 0; i < free_rates.size(); ++i) {
    sum += free_rates[i] * _rate_weights[p_index][i];
  }
  for (auto &f : free_rates) { f /= sum; }
  debug_print(
      EMIT_LEVEL_DEBUG, "setting rates to %s", to_string(free_rates).c_str());
  pll_set_category_rates(_partitions[p_index], _rate_rates[p_index].data());
}

void model_t::update_invariant_sites(size_t p_index) {
  if (_invariant_sites) {
    pll_update_invariant_sites(_partitions[p_index]);
  } else {
    for (unsigned int i = 0; i < _submodels; ++i) {
      pll_update_invariant_sites_proportion(_partitions[p_index], i, 0.0);
    }
  }
}

void model_t::set_tip_states(size_t p_index, const msa_t &msa) {
  /* make a label map */
  auto label_map = _tree.label_map();

  /* use the label map to assign tip states in the partition */

  for (int i = 0; i < msa.count(); ++i) {
    try {
      auto result = pll_set_tip_states(_partitions[p_index],
                                       label_map.at(msa.label(i)),
                                       msa.map(),
                                       msa.sequence(i));
      if (result == PLL_FAILURE) {
        throw std::runtime_error("failed to set tip " + std::to_string(i));
      }
    } catch (const std::exception &e) {
      throw std::runtime_error(std::string("Could not find taxa ")
                               + msa.label(i) + " in tree");
    }
  }

  /* set pattern weights */
  pll_set_pattern_weights(_partitions[p_index], msa.weights());
}

void model_t::set_empirical_freqs(size_t p_index) {
  pll_partition_t *partition = _partitions[p_index];
  double *         emp_freqs = pllmod_msa_empirical_frequencies(partition);
  for (size_t i = 0; i < partition->states; ++i) {
    if (emp_freqs[i] <= 0) {
      throw invalid_empirical_frequencies_exception(
          "One of the state frequenices is zero while using emperical "
          "frequencies");
    }
  }
  pll_set_frequencies(partition, 0, emp_freqs);
  free(emp_freqs);
}

void model_t::set_freqs(size_t p_index, const model_params_t &freqs) {
  for (auto f : freqs) {
    if (f <= 0.0) {
      throw std::runtime_error("Frequencies with 0 entries are not allowed");
    }
  }
  pll_set_frequencies(_partitions[p_index], 0, freqs.data());
}

void model_t::set_freqs_all_free(size_t p_index, model_params_t freqs) {
  double sum = 0.0;
  for (auto p : freqs) { sum += p; }
  for (auto &f : freqs) { f /= sum; }
  set_freqs(p_index, freqs);
}

bool model_t::update_eigen_partition(size_t partition_index) {
  bool updated = false;
  auto part    = _partitions[partition_index];
  for (size_t i = 0; i < part->rate_cats; ++i) {
    if (!part->eigen_decomp_valid[i]) {
      int result = pll_update_eigen(part, _param_indicies[partition_index][i]);
      updated    = true;
      if (result == PLL_FAILURE) {
        throw std::runtime_error{"Failed to update the eigen decomposition"};
      }
    }
  }
  return updated;
}

void model_t::update_pmatrix_partition(
    size_t                           partition_index,
    const std::vector<unsigned int> &pmatrix_indices,
    const std::vector<double> &      branch_lengths) {
  auto part = _partitions[partition_index];
#pragma omp parallel for collapse(1) schedule(static)
  for (size_t branch = 0; branch < branch_lengths.size(); ++branch) {
    auto param_index   = _param_indicies[partition_index];
    auto matrix_index  = pmatrix_indices[branch];
    auto branch_length = branch_lengths[branch];
    pll_update_prob_matrices(
        part, param_index.data(), &matrix_index, &branch_length, 1);
  }
}

std::vector<bool>
model_t::update_pmatrices(const std::vector<unsigned int> &pmatrix_indices,
                          const std::vector<double> &      branch_lengths) {
  /* Update the eigen decompositions first */
  std::vector<bool> updated(_partitions.size(), false);
  for (size_t part_index = 0; part_index < _partitions.size(); ++part_index) {
    updated[part_index] = update_eigen_partition(part_index);
  }

  for (size_t part_index = 0; part_index < _partitions.size(); ++part_index) {
    update_pmatrix_partition(part_index, pmatrix_indices, branch_lengths);
  }
  return updated;
}

double model_t::compute_lh(const root_location_t &root_location) {
  std::vector<pll_operation_t> ops;
  std::vector<unsigned int>    pmatrix_indices;
  std::vector<double>          branch_lengths;
  bool new_root = root_location != _tree.root_location();

  GENERATE_AND_UNPACK_OPS(
      _tree, root_location, ops, pmatrix_indices, branch_lengths);

  auto updated_partitions = update_pmatrices(pmatrix_indices, branch_lengths);

  double lh = 0.0;

#pragma omp parallel for reduction(+ : lh)
  for (size_t i = 0; i < _partitions.size(); ++i) {
    auto &partition = _partitions[i];

    if (new_root || updated_partitions[i]) {
      pll_update_partials(
          partition, ops.data(), static_cast<unsigned int>(ops.size()));
    }

    lh += pll_compute_root_loglikelihood(partition,
                                         _tree.root_clv_index(),
                                         _tree.root_scaler_index(),
                                         _param_indicies[i].data(),
                                         nullptr);
  }
  return lh;
}

double model_t::compute_lh_root(const root_location_t &root) {
  pll_operation_t           op;
  std::vector<unsigned int> matrix_indices;
  std::vector<double>       branch_lengths;

  {
    auto result    = _tree.generate_derivative_operations(root);
    op             = std::move(std::get<0>(result));
    matrix_indices = std::move(std::get<1>(result));
    branch_lengths = std::move(std::get<2>(result));
  }

  double lh = 0.0;

#pragma omp parallel for reduction(+ : lh)
  for (size_t i = 0; i < _partitions.size(); ++i) {
    auto &partition = _partitions[i];
    int   result    = pll_update_prob_matrices(
        partition,
        _param_indicies[i].data(),
        matrix_indices.data(),
        branch_lengths.data(),
        static_cast<unsigned int>(matrix_indices.size()));

    if (result == PLL_FAILURE) { throw std::runtime_error(pll_errmsg); }
    pll_update_partials(partition, &op, 1);
    lh += pll_compute_root_loglikelihood(partition,
                                         _tree.root_clv_index(),
                                         _tree.root_scaler_index(),
                                         _param_indicies[i].data(),
                                         nullptr);
    //* _partition_weights[i];
  }
  if (std::isnan(lh)) {
    throw std::runtime_error("lh at root is not a number: "
                             + std::to_string(lh));
  }
  return lh;
}

double
model_t::compute_lh_partition(size_t partition_index,
                              const std::vector<pll_operation_t> &ops,
                              const std::vector<unsigned int> &pmatrix_indices,
                              const std::vector<double> &      branch_lengths) {
  bool update_partials = update_eigen_partition(partition_index);
  if (update_partials) {
    update_pmatrix_partition(partition_index, pmatrix_indices, branch_lengths);
    pll_update_partials(_partitions[partition_index],
                        ops.data(),
                        static_cast<unsigned int>(ops.size()));
  }

  double lh =
      pll_compute_root_loglikelihood(_partitions[partition_index],
                                     _tree.root_clv_index(),
                                     _tree.root_scaler_index(),
                                     _param_indicies[partition_index].data(),
                                     nullptr);
  if (std::isnan(lh)) {
    throw std::runtime_error("lh at root is not a number: "
                             + std::to_string(lh));
  }
  return lh;
}

/*
 * Use a secant method to compute the derivative
 */
dlh_t model_t::compute_dlh(const root_location_t &root) {
  dlh_t            ret;
  constexpr double EPSILON = 1e-8;
  root_location_t  root_prime{root};
  root_prime.brlen_ratio += EPSILON;
  double sign = 1.0;
  if (root_prime.brlen_ratio >= 1.0) {
    root_prime.brlen_ratio = root.brlen_ratio - EPSILON;
    sign                   = -1.0;
  }

  double fx = compute_lh_root(root);
  ret.lh    = fx;

  if (std::isnan(fx)) {
    throw std::runtime_error("fx is not finite when computing derivative: "
                             + std::to_string(root.edge->length));
  }
  double fxh = compute_lh_root(root_prime);
  if (std::isnan(fxh)) {
    throw std::runtime_error("fxh is not finite when computing derivative: "
                             + std::to_string(root_prime.edge->length));
  }

  if (std::isnan(fxh)) {
    throw std::runtime_error("fxh is not finite when computing derivative: "
                             + std::to_string(root_prime.edge->length));
  }
  if (std::isinf(fxh) && std::isinf(fx)) {
    debug_string(EMIT_LEVEL_DEBUG,
                 "Both evals are -inf, returning a 0 derivative");
    return {fx, 0};
  }
  double dlh = (fxh - fx) / EPSILON;
  debug_print(
      EMIT_LEVEL_DEBUG, "dlh: %f, fx: %f, fxh: %f", dlh * sign, fx, fxh);
  ret.dlh = dlh * sign;
  return ret;
}

/*
 * Find the root of dlh w.r.t. alpha in order to find the optimum value.
 * Technically, this also evalutes lh to find the true maximum, so it isn't
 * strictly a bisection method.
 *
 * Some tunables:
 *  - depth
 *
 * Generally the prefered function for this is brents. This is left here for
 * legacy reasons.
 */
std::pair<root_location_t, double> model_t::bisect(const root_location_t &beg,
                                                   dlh_t                  d_beg,
                                                   const root_location_t &end,
                                                   dlh_t                  d_end,
                                                   double                 atol,
                                                   size_t depth = 0) {
  assert_string(d_beg.dlh * d_end.dlh > 0,
                "Bisect called with endpoints which don't bracket");
  root_location_t midpoint{beg};
  midpoint.brlen_ratio = (beg.brlen_ratio + end.brlen_ratio) / 2;

  auto d_midpoint = compute_dlh(midpoint);

  if (depth > 32) {
    debug_print(
        EMIT_LEVEL_DEBUG,
        "depth exceeded limit, returning midpoint with ratio: %f, lh: %f",
        midpoint.brlen_ratio,
        d_midpoint.lh);
    return {midpoint, d_midpoint.lh};
  }

  /* case 1: d_midpoint is within tolerance of 0.0
   * We found a root, and should return
   */

  if (fabs(d_midpoint.dlh) < atol) {
    debug_string(EMIT_LEVEL_DEBUG, "case 1");
    return {midpoint, d_midpoint.lh};
  }

  /* case 2: midpoint is opposite sign of both beg and end
   * There are at least 2 roots, so we recurse on both sides, and return the
   * one with the best LH
   */
  if ((d_beg.dlh > 0.0 && d_end.dlh > 0.0 && d_midpoint.dlh < 0.0)
      || (d_beg.dlh < 0.0 && d_end.dlh < 0.0 && d_midpoint.dlh > 0.0)) {
    debug_string(EMIT_LEVEL_DEBUG, "case 2");
    auto r1 = bisect(beg, d_beg, midpoint, d_midpoint, atol, depth + 1);
    auto r2 = bisect(midpoint, d_midpoint, end, d_end, atol, depth + 1);
    if (r1.second < r2.second) { return r2; }
    return r1;
  }
  /* case 3: end and midpoint share a sign, while beg has the opposite sign
   * In this case, there is a root between midpoint and beg
   */
  if ((d_beg.dlh < 0.0 && d_midpoint.dlh > 0.0 && d_end.dlh > 0.0)
      || (d_beg.dlh > 0.0 && d_midpoint.dlh < 0.0 && d_end.dlh < 0.0)) {
    debug_print(
        EMIT_LEVEL_DEBUG,
        "case 3, d_beg.dlh: %f, d_midpoint.dlh: %f, d_end:%f, alpha: %f",
        d_beg.lh,
        d_midpoint.dlh,
        d_end.dlh,
        midpoint.brlen_ratio);
    return bisect(beg, d_beg, midpoint, d_midpoint, atol, depth + 1);
  }

  /* case 4: beg and midpoint share a sign, while end has the opposite sign
   * In this case, there is a root between midpoint and end
   */
  if ((d_beg.dlh < 0.0 && d_midpoint.dlh < 0.0 && d_end.dlh > 0.0)
      || (d_beg.dlh > 0.0 && d_midpoint.dlh > 0.0 && d_end.dlh < 0.0)) {
    debug_string(EMIT_LEVEL_DEBUG, "case 4");
    return bisect(midpoint, d_midpoint, end, d_end, atol, depth + 1);
  }

  /* case 5: something went wrong */
  throw std::runtime_error("Bisection failed to converge with interval : "
                           + std::to_string(d_beg.dlh) + ", "
                           + std::to_string(d_end.dlh) + ", midpoint.dlh: "
                           + std::to_string(d_midpoint.dlh));
}

std::pair<root_location_t, double> model_t::brents(root_location_t beg,
                                                   dlh_t           d_beg,
                                                   root_location_t end,
                                                   dlh_t           d_end,
                                                   double          atol) {
  assert_string(d_beg.dlh * d_end.dlh < 0,
                "Brents called with endpoints which don't bracket");

  root_location_t midpoint{end};
  auto            d_midpoint = d_end;
  double          e;
  double          d;
  d = e = end.brlen_ratio - beg.brlen_ratio;

  for (size_t i = 0; i < 64; ++i) {
    if (d_end.dlh * d_midpoint.dlh > 0.0) {
      midpoint   = beg;
      d_midpoint = d_beg;
      d = e = end.brlen_ratio - beg.brlen_ratio;
    }
    if (fabs(d_end.dlh) < fabs(d_midpoint.dlh)) {
      beg        = end;
      end        = midpoint;
      midpoint   = beg;
      d_beg      = d_end;
      d_end      = d_midpoint;
      d_midpoint = d_beg;
    }

    double tol =
        2.0 * fabs(end.brlen_ratio) * std::numeric_limits<double>::epsilon()
        + 0.5 * atol;
    double e_tol = 0.5 * (midpoint.brlen_ratio - end.brlen_ratio);
    if (fabs(e_tol) <= tol || fabs(d_end.dlh) <= 1e-12) return {end, d_end.lh};
    if (fabs(e) >= tol && fabs(d_beg.dlh) > fabs(d_end.dlh)) {
      double s = d_end.dlh / d_beg.dlh;
      double p, q;
      if (fabs(beg.brlen_ratio - midpoint.brlen_ratio) < 1e-12) {
        p = 2.0 * e_tol * s;
        q = 1.0 - s;
      } else {
        q        = d_beg.dlh / d_midpoint.dlh;
        double r = d_end.dlh / d_midpoint.dlh;
        p        = s
            * (2.0 * e_tol * q * (q - r)
               - (end.brlen_ratio - beg.brlen_ratio) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      p           = fabs(p);
      double min1 = 3.0 * e_tol * q - fabs(e_tol * q);
      double min2 = fabs(e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        e = d;
        d = p / q;
      } else {
        d = e_tol;
        e = d;
      }
    } else {
      d = e_tol;
      e = d;
    }
    beg   = end;
    d_beg = d_end;
    if (fabs(d) > tol) end.brlen_ratio += d;
    else {
      end.brlen_ratio += e_tol >= 0.0 ? tol : -tol;
    }
    d_end = compute_dlh(end);
  }
  throw std::runtime_error("Brents method failed to converge");
}

/* Find the optimum for the ratio via brents method */
root_location_t model_t::optimize_alpha(const root_location_t &root,
                                        double                 atol) {
  double lh = compute_lh_root(root);
  if (std::isnan(lh)) {
    throw std::runtime_error("initial likelihood calculation is not finite");
  }
  root_location_t beg{root};
  beg.brlen_ratio = 0.0;

  root_location_t end{root};
  end.brlen_ratio = 1.0;

  auto d_beg = compute_dlh(beg);

  auto d_end = compute_dlh(end);

  if (std::isnan(d_beg.dlh) || std::isnan(d_end.dlh)) {
    throw std::runtime_error(
        "Initial derivatives failed when optimizing alpha: "
        + std::to_string(root.edge->length));
  }

  root_location_t best_endpoint    = d_beg.lh >= d_end.lh ? beg : end;
  auto            lh_best_endpoint = d_beg.lh >= d_end.lh ? d_beg : d_end;
  if (d_beg.lh >= d_end.lh) {
    debug_string(EMIT_LEVEL_DEBUG, "beg endpoint is best");
  } else {
    debug_string(EMIT_LEVEL_DEBUG, "end endpoint is best");
  }

  debug_print(EMIT_LEVEL_DEBUG, "lh_endpoint.dlh: %f", lh_best_endpoint.dlh);
  debug_print(
      EMIT_LEVEL_DEBUG, "d_beg.dlh: %f, d_end.dlh: %f", d_beg.dlh, d_end.dlh);

  if (fabs(d_beg.dlh) < atol || fabs(d_end.dlh) < atol) {
    debug_string(EMIT_LEVEL_DEBUG, "one of the endpoints is sufficient");
    return best_endpoint;
  }

  if ((d_beg.dlh < 0.0 && d_end.dlh > 0.0)
      || (d_beg.dlh > 0.0 && d_end.dlh < 0.0)) {
    auto mid = brents(beg, d_beg, end, d_end, atol);
    debug_print(EMIT_LEVEL_DEBUG,
                "mid lh: %f, end lh: %f",
                mid.second,
                lh_best_endpoint.lh);
    return lh_best_endpoint.lh > mid.second ? best_endpoint : mid.first;
  }

  /* if we have gotten here, then we must have an "even" function so, grid
   * search to find a midpoint which has an opposite sign value and use that
   * to do bisection with
   */

  bool  beg_end_pos      = d_beg.dlh > 0.0 && d_end.dlh > 0.0;
  dlh_t best_midpoint_lh = {-std::numeric_limits<double>::infinity(), 0};
  root_location_t best_midpoint;
  bool            found_midpoint = false;

  for (size_t midpoints = 2; midpoints <= 32; midpoints *= 2) {
    for (size_t midpoint = 1; midpoint <= midpoints; ++midpoint) {
      if (midpoint % 2 == 0) continue;

      double alpha = 1.0 / (double)midpoints * midpoint;
      debug_print(EMIT_LEVEL_DEBUG, "alpha: %f", alpha);
      root_location_t midpoint_root{beg};
      midpoint_root.brlen_ratio = alpha;
      auto d_midpoint           = compute_dlh(midpoint_root);
      debug_print(EMIT_LEVEL_DEBUG, "d_midpoint.dlh: %f", d_midpoint.dlh);
      if (fabs(d_midpoint.dlh) < atol) {
        if (best_midpoint_lh.lh < d_midpoint.lh) {
          best_midpoint_lh = d_midpoint;
          best_midpoint    = midpoint_root;
          found_midpoint   = true;
        }
      }
      if ((beg_end_pos && d_midpoint.dlh < 0.0)
          || (!beg_end_pos && d_midpoint.dlh > 0.0)) {
        /*
         * we have a midpoint, so now we need to figure out if it is a min or
         * a maximum.
         */
        auto r1 = brents(beg, d_beg, midpoint_root, d_midpoint, atol);
        auto r2 = brents(midpoint_root, d_midpoint, end, d_end, atol);
        debug_print(
            EMIT_LEVEL_DEBUG, "r1 lh: %f, r2 lh: %f", r1.second, r2.second);
        if (lh_best_endpoint.lh < best_midpoint_lh.lh) {
          lh_best_endpoint = best_midpoint_lh;
          best_endpoint    = best_midpoint;
        }
        if (r1.second < r2.second) {
          return lh_best_endpoint.lh >= r2.second ? best_endpoint : r2.first;
        }
        return lh_best_endpoint.lh >= r1.second ? best_endpoint : r1.first;
      }
    }
  }

  if (found_midpoint) { return best_midpoint; }

  /*
   * If we got here, we can just return the best of beg or end, since it seems
   * like there is no root. If both beg and end are positive, then that means
   * the liklihood _increases_ as we increase alpha, so the best liklihood is
   * at the end. Likewise, if both are negative, then the best endpoint is at
   * the begining.
   */
  if (beg_end_pos) {
    return end;
  } else {
    return beg;
  }

  throw std::runtime_error(
      "Initial derivatives failed when optimizing alpha, ran out of cases");
}

std::pair<root_location_t, double>
model_t::optimize_root_location(size_t min_roots, double root_ratio) {
  std::pair<root_location_t, double> best;
  best.second = -std::numeric_limits<double>::infinity();

  /* start by making a list of "good" roots, with the current model*/

  auto sorted_roots = suggest_roots_random(min_roots, root_ratio);
  for (auto &sr : sorted_roots) {
    auto &rl = sr.first;
    debug_print(EMIT_LEVEL_DEBUG, "working rl: %s", rl.label().c_str());

    move_root(rl);
    rl = optimize_alpha(rl, 1e-14);
    debug_print(EMIT_LEVEL_DEBUG, "alpha: %f", rl.brlen_ratio);

    double rl_lh = compute_lh_root(rl);
    debug_print(EMIT_LEVEL_DEBUG, "rl_lh: %f", rl_lh);

    if (rl_lh > best.second) {
      best.first  = rl;
      best.second = rl_lh;
    }
  }
  debug_print(EMIT_LEVEL_DEBUG, "finished with lh: %f", best.second);
  return best;
}

void model_t::move_root(const root_location_t &new_root) {
  std::vector<pll_operation_t> ops;
  std::vector<unsigned int>    pmatrix_indices;
  std::vector<double>          branch_lengths;

  {
    auto results    = _tree.generate_root_update_operations(new_root);
    ops             = std::move(std::get<0>(results));
    pmatrix_indices = std::move(std::get<1>(results));
    branch_lengths  = std::move(std::get<2>(results));
  }
  debug_print(EMIT_LEVEL_DEBUG,
              "ops.size: %lu, pmats.size: %lu, brlens.size: %lu",
              ops.size(),
              pmatrix_indices.size(),
              branch_lengths.size());

  for (size_t i = 0; i < _partitions.size(); ++i) {
    auto &partition = _partitions[i];
    int   result    = pll_update_prob_matrices(
        partition,
        _param_indicies[i].data(),
        pmatrix_indices.data(),
        branch_lengths.data(),
        static_cast<unsigned int>(pmatrix_indices.size()));

    if (result == PLL_FAILURE) { throw std::runtime_error(pll_errmsg); }

    pll_update_partials(
        partition, ops.data(), static_cast<unsigned int>(ops.size()));
  }
}

std::vector<std::pair<root_location_t, double>> model_t::suggest_roots() {
  return suggest_roots_random(1, 0.0);
}

std::vector<std::pair<root_location_t, double>>
model_t::suggest_roots_random(size_t min, double ratio) {
  std::vector<std::pair<root_location_t, double>> rl_lhs;
  using fucking_difference_type =
      std::vector<std::pair<root_location_t, double>>::difference_type;
  rl_lhs.reserve(_tree.root_count());
  for (auto rl : _tree.roots()) {
    move_root(rl);
    rl_lhs.push_back(std::make_pair(rl, compute_lh_root(rl)));
  }
  auto final_size = std::max(static_cast<size_t>(rl_lhs.size() * ratio), min);
  std::partial_sort(
      rl_lhs.begin(),
      rl_lhs.begin() + static_cast<fucking_difference_type>(final_size),
      rl_lhs.end(),
      [](std::pair<root_location_t, double> a,
         std::pair<root_location_t, double> b) { return a.second > b.second; });
  rl_lhs.resize(final_size);
  return rl_lhs;
}

std::vector<std::pair<root_location_t, double>>
model_t::suggest_roots_midpoint(size_t min, double ratio) {
  auto midpoints = _tree.rank_midpoints();
  auto final_size =
      std::max(static_cast<size_t>(midpoints.size() * ratio), min);
  std::vector<std::pair<root_location_t, double>> ret;
  ret.reserve(final_size);
  for (size_t i = 0; i < final_size; i++) {
    ret.push_back(std::make_pair(midpoints[i], compute_lh(midpoints[i])));
  }
  return ret;
}

std::vector<size_t> model_t::shuffle_root_indicies() {
  std::vector<size_t> idx(_tree.root_count());
  std::iota(idx.begin(), idx.end(), 0);
  std::shuffle(idx.begin(), idx.end(), _random_engine);
  return idx;
}

partition_parameters_t model_t::make_partition_parameters(
    size_t states, rate_category::rate_category_e rc, size_t rate_cat_count) {
  partition_parameters_t pp;
  size_t                 subst_size = states * states - states;

  pp.subst_rates.resize(subst_size);
  for (auto &v : pp.subst_rates) { v = 1.0 / subst_size; }

  pp.freqs.resize(states);
  std::fill(pp.freqs.begin(), pp.freqs.end(), 1.0 / states);
  for (auto &v : pp.freqs) { v = 1.0 / states; }

  switch (rc) {
  case rate_category::MEAN:
  case rate_category::MEDIAN:
    pp.gamma_alpha.resize(1);
    pp.gamma_alpha[0] = 1.0;
    break;
  case rate_category::FREE:
    pp.gamma_alpha.resize(rate_cat_count);
    for (auto &v : pp.gamma_alpha) { v = 1.0; }
    pp.gamma_weights.resize(rate_cat_count);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (auto &v : pp.gamma_weights) { v = dis(_random_engine); }
  }
  return pp;
}

/* Optimize the substitution parameters and root location.*/
std::pair<root_location_t, double> model_t::search(size_t        min_roots,
                                                   double        root_ratio,
                                                   double        atol,
                                                   double        pgtol,
                                                   double        brtol,
                                                   double        factor,
                                                   checkpoint_t &checkpoint) {
  if (_assigned_idx.size() == 0) {
    debug_string(EMIT_LEVEL_WARNING, "There is no work to be done");
  }
  double          best_lh = -std::numeric_limits<double>::infinity();
  root_location_t best_rl;

  set_subst_rates_uniform();
  set_empirical_freqs();
  size_t root_count = _assigned_idx.size();
  size_t root_index = 0;

  std::vector<partition_parameters_t> best_params;
  best_params.reserve(_partitions.size());

  debug_string(EMIT_LEVEL_PROGRESS, "Starting root search");

  for (auto rl_index : _assigned_idx) {
    auto rl = _tree.root_location(rl_index);
    set_subst_rates_uniform();
    set_empirical_freqs();

    std::vector<partition_parameters_t> params;
    params.reserve(_partitions.size());

    std::vector<partition_parameters_t> saved_params;
    saved_params.reserve(_partitions.size());

    for (size_t p = 0; p < _partitions.size(); ++p) {
      params.push_back(make_partition_parameters(_partitions[p]->states,
                                                 _rate_category_types[p],
                                                 _partitions[p]->rate_cats));
    }

    auto   cur_best_rl = rl;
    double cur_best_lh = -std::numeric_limits<double>::infinity();

    for (size_t iter = 0; iter < 1e3; ++iter) {
      saved_params = params;

      optimize_params(params, rl, pgtol, factor, true);

      debug_string(EMIT_LEVEL_INFO, "Optimizing Root Location");
      auto cur = optimize_root_location(min_roots, root_ratio);

      debug_print(EMIT_LEVEL_INFO, "Iteration %lu LH: %.9f", iter, cur.second);

      if (cur.second < cur_best_lh) {
        /* We failed to make any progress, so just give up */
        debug_string(EMIT_LEVEL_DEBUG,
                     "breaking due to failure to make progress");
        for (size_t i = 0; i < _partitions.size(); ++i) {
          set_subst_rates(i, saved_params[i].subst_rates);
          set_freqs(i, saved_params[i].freqs);
          set_gamma_rates(i, saved_params[i].gamma_alpha);
          if (_rate_category_types[i] == rate_category::FREE) {
            set_gamma_weights(i, saved_params[i].gamma_weights);
          }
        }
        params = saved_params;
        break;
      }

      if (_early_stop) {
        if (rl.edge == cur.first.edge
            && fabs(rl.brlen_ratio - cur.first.brlen_ratio) < brtol) {
          debug_string(EMIT_LEVEL_DEBUG, "breaking due to early stop");
          cur_best_rl = cur.first;
          cur_best_lh = cur.second;
          break;
        }
      }

      if (fabs(cur.second - cur_best_lh) < atol) {
        debug_string(EMIT_LEVEL_DEBUG, "breaking due to atol");
        cur_best_rl = cur.first;
        cur_best_lh = cur.second;
        break;
      }

      cur_best_rl = cur.first;
      cur_best_lh = cur.second;

      rl = cur_best_rl;
    }

    ++root_index;
    debug_print(EMIT_LEVEL_PROGRESS,
                "Stage %lu/%lu, ETC: %fh",
                root_index,
                root_count,
                progress_macro(root_index, root_count));

    checkpoint.write({cur_best_rl.id, cur_best_lh, cur_best_rl.brlen_ratio},
                     params);

    debug_print(EMIT_LEVEL_DEBUG,
                "finished optimize_all root, cur_best_lh: %f",
                cur_best_lh);
  }

#ifdef MPI_VERSION
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (__MPI_RANK__ == 0) {
    auto total_progress    = checkpoint.read_results();
    auto total_best_result = *std::max_element(
        total_progress.begin(),
        total_progress.end(),
        [](std::pair<rd_result_t, std::vector<partition_parameters_t>> a,
           std::pair<rd_result_t, std::vector<partition_parameters_t>> b)
            -> bool { return a.first.lh < b.first.lh; });

    best_rl             = _tree.root_location(total_best_result.first.root_id);
    best_rl.brlen_ratio = total_best_result.first.alpha;
    best_lh             = total_best_result.first.lh;
    best_params         = total_best_result.second;
    set_model_params(best_params);
  }

  move_root(best_rl);
  return {best_rl, best_lh};
}

std::pair<root_location_t, double>
model_t::exhaustive_search(double        atol,
                           double        pgtol,
                           double        brtol,
                           double        factor,
                           checkpoint_t &checkpoint) {
  if (_assigned_idx.size() == 0) {
    debug_string(EMIT_LEVEL_WARNING, "There is no work to be done");
  }
  size_t          root_index = 0;
  root_location_t best_rl;
  double          best_lh    = -std::numeric_limits<double>::infinity();
  size_t          root_count = _assigned_idx.size();
  debug_string(EMIT_LEVEL_PROGRESS, "Starting exhaustive search");

  for (auto rl_index : _assigned_idx) {
    auto rl = _tree.root_location(rl_index);
    set_subst_rates_uniform();
    set_empirical_freqs();

    _tree.root_by(rl);
    compute_lh(rl);
    std::vector<partition_parameters_t> params;

    for (size_t p = 0; p < _partitions.size(); ++p) {
      params.push_back(make_partition_parameters(_partitions[p]->states,
                                                 _rate_category_types[p],
                                                 _partitions[p]->rate_cats));
    }

    root_location_t cur_best_rl = rl;
    double          cur_best_lh = -std::numeric_limits<double>::infinity();

    for (size_t iter = 0; iter < 1e3; ++iter) {
      debug_string(EMIT_LEVEL_MPROGRESS, "Optimizing parameters");
      optimize_params(params, rl, pgtol, factor, (iter % 10 == 0));

      if (fabs(compute_lh(rl) - cur_best_lh) < atol) { break; }

      debug_string(EMIT_LEVEL_MPROGRESS, "Optimizing Root Location");
      auto   cur_rl = optimize_alpha(rl, brtol);
      double cur_lh = compute_lh_root(cur_rl);

      debug_print(EMIT_LEVEL_MPROGRESS, "Iteration %lu LH: %.5f", iter, cur_lh);
      debug_print(
          EMIT_LEVEL_INFO, "difference in lh: %.5f", (cur_lh - cur_best_lh));

      if (_early_stop) {
        if (fabs(rl.brlen_ratio - cur_rl.brlen_ratio) < brtol) {
          debug_print(EMIT_LEVEL_DEBUG,
                      "Current BRlen ratio tolerances: %.7f, brtol: %.7f",
                      fabs(rl.brlen_ratio - cur_rl.brlen_ratio),
                      brtol);
          cur_best_rl = cur_rl;
          cur_best_lh = cur_lh;
          break;
        }
      }

      if ((cur_lh - cur_best_lh) < atol) {
        if (cur_lh > cur_best_lh) {
          cur_best_rl = cur_rl;
          cur_best_lh = cur_lh;
        }
        break;
      }

      if (cur_lh > cur_best_lh) {
        cur_best_rl = cur_rl;
        cur_best_lh = cur_lh;
      }

      rl = cur_rl;
    }

    checkpoint.write({cur_best_rl.id, cur_best_lh, cur_best_rl.brlen_ratio},
                     params);
    root_index++;

    debug_print(EMIT_LEVEL_PROGRESS,
                "Step %lu / %lu, ETC: %0.2fh",
                root_index,
                root_count,
                progress_macro(root_index, root_count));

    if (cur_best_lh > best_lh) {
      best_rl = cur_best_rl;
      best_lh = cur_best_lh;
    }
  }

#ifdef MPI_VERSION
  debug_string(EMIT_LEVEL_IMPORTANT, "Waiting for the rest to finish");
  MPI_Barrier(MPI_COMM_WORLD);
  debug_string(EMIT_LEVEL_IMPORTANT, "Done waiting");
#endif

  if (__MPI_RANK__ == 0) {
    auto   total_progress = checkpoint.read_results();
    double max_lh         = -std::numeric_limits<double>::infinity();

    for (auto result : total_progress) {
      max_lh = std::max(result.first.lh, max_lh);
    }

    double total_lh = 0;
    for (auto result : total_progress) {
      total_lh += exp(result.first.lh - max_lh);
    }
    debug_print(EMIT_LEVEL_DEBUG, "LWR denom: %f, %e", total_lh, total_lh - 1);

    for (auto result : total_progress) {
      double lwr     = exp((result.first.lh - max_lh)) / total_lh;
      auto   rl      = _tree.root_location(result.first.root_id);
      rl.brlen_ratio = result.first.alpha;
      _tree.annotate_branch(rl, "LWR", std::to_string(lwr));
      _tree.annotate_lh(rl, result.first.lh);
      _tree.annotate_ratio(rl, result.first.alpha);
    }
    auto total_best_result = *std::max_element(
        total_progress.begin(),
        total_progress.end(),
        [](std::pair<rd_result_t, std::vector<partition_parameters_t>> a,
           std::pair<rd_result_t, std::vector<partition_parameters_t>> b)
            -> bool { return a.first.lh < b.first.lh; });

    best_rl             = _tree.root_location(total_best_result.first.root_id);
    best_rl.brlen_ratio = total_best_result.first.alpha;
    best_lh             = total_best_result.first.lh;
  }

  return {best_rl, best_lh};
}

void model_t::initialize() { compute_lh(_tree.root_location(0)); }

void model_t::finalize() { _tree.unroot(); }

rooted_tree_t model_t::rooted_tree(const root_location_t &root) const {
  rooted_tree_t new_tree(_tree);
  new_tree.root_by(root.id);
  return new_tree;
}

rooted_tree_t model_t::virtual_rooted_tree(const root_location_t &root) const {
  rooted_tree_t new_tree(_tree);
  new_tree.root_by(root.id);
  new_tree.unroot();
  return new_tree;
}

rooted_tree_t model_t::unrooted_tree() const {
  rooted_tree_t new_tree(_tree);
  new_tree.unroot();
  return new_tree;
}

void model_t::initialize_partitions(const std::vector<msa_t> &msa) {
  for (size_t partition_index = 0; partition_index < _partitions.size();
       ++partition_index) {
    set_tip_states(partition_index, msa[partition_index]);
    update_invariant_sites(partition_index);
    set_empirical_freqs(partition_index);
    set_subst_rates_random(partition_index, msa[partition_index]);
    // set_gamma_rates(partition_index);
  }
}

void model_t::initialize_partitions_uniform_freqs(
    const std::vector<msa_t> &msa) {
  for (size_t partition_index = 0; partition_index < _partitions.size();
       ++partition_index) {
    set_tip_states(partition_index, msa[partition_index]);
    update_invariant_sites(partition_index);
    std::vector<double> uni_freqs(
        _partitions[partition_index]->states,
        1.0 / (double)_partitions[partition_index]->states);
    set_freqs(partition_index, uni_freqs);
    set_subst_rates_random(partition_index, msa[partition_index]);
    set_gamma_rates(partition_index);
  }
}

std::string model_t::subst_string() const {
  std::ostringstream oss;
  oss << "{";
  for (size_t p_index = 0; p_index < _partitions.size(); ++p_index) {
    auto p = _partitions[p_index];
    oss << "{";
    size_t              params_size = p->states * p->states - p->states;
    std::vector<double> params{p->subst_params[0],
                               p->subst_params[0] + params_size};
    for (size_t i = 0; i < params.size(); ++i) {
      oss << std::to_string(params[i]);
      if (i != params.size() - 1) oss << ",";
    }
    oss << "}";
    if (p_index != _partitions.size() - 1) oss << ",";
  }
  oss << "}";
  return oss.str();
}

static double
gd_params(model_params_t &                                    initial_params,
          size_t                                              partition_index,
          double                                              p_min,
          double                                              p_max,
          double                                              epsilon,
          std::function<double()>                             compute_lh,
          std::function<void(size_t, const model_params_t &)> set_func) {
  size_t              iters    = 0;
  double              atol     = 1e-4;
  size_t              n_params = initial_params.size();
  std::vector<double> parameters(initial_params);
  std::vector<double> gradient(n_params, 0.0);
  set_func(partition_index, initial_params);
  double initial_score = compute_lh();
  double last_score    = std::numeric_limits<double>::infinity();
  while (true) {
    set_func(partition_index, parameters);
    double score = compute_lh();
    debug_print(EMIT_LEVEL_INFO, "GD Iter: %lu Score: %.5f", iters, -score);
    if (fabs(last_score - score) < atol) { break; }
    for (size_t i = 0; i < n_params; ++i) {
      double h = epsilon * fabs(parameters[i]);
      if (h < epsilon) { h = epsilon; }
      double temp = parameters[i];
      parameters[i] += h;
      set_func(partition_index, parameters);
      double dlh = compute_lh();
      assert_string(std::isfinite(dlh), "dlh is not finite");
      gradient[i] = (dlh - score) / h;
      assert_string(std::isfinite(gradient[i]), "gradient is not finite");
      parameters[i] = temp;
    }
    /* backtracking line search */
    double gnorm = 0.0;
    /* compute the norm of the gradient */
    for (size_t j = 0; j < n_params; ++j) {
      gnorm += gradient[j] * gradient[j];
    }
    gnorm = sqrt(gnorm);
    /* compute p * gnorm */
    double pgnorm = 0.0;
    for (size_t j = 0; j < n_params; ++j) {
      pgnorm += gradient[j] * gradient[j] / gnorm;
    }
    /* search for the correct alpha */
    double alpha = 1.0;
    while (alpha > 1e-12) {
      double q = score + 1e-4 * alpha * pgnorm;
      auto   tmp_params{parameters};
      for (size_t j = 0; j < n_params; ++j) {
        tmp_params[j] = parameters[j] - alpha * gradient[j];
        if (tmp_params[j] < p_min) {
          tmp_params[j] = p_min;
        } else if (tmp_params[j] > p_max) {
          tmp_params[j] = p_max;
        }
      }
      set_func(partition_index, tmp_params);
      double s = compute_lh();
      if (q > s) {
        std::swap(tmp_params, parameters);
        break;
      }
      alpha /= 2;
    }

    last_score = score;
    set_func(partition_index, parameters);
    iters++;
  }
  set_func(partition_index, parameters);
  double final_score = compute_lh();
  // check that we have successfully _minimized_ the function
  // assert_string(initial_score >= final_score, "Failed to make the lh
  // better");
  if (initial_score >= final_score) {
    std::swap(parameters, initial_params);
  } else {
    debug_print(EMIT_LEVEL_WARNING,
                "Failed to improve the likelihood after GD. Start: %f, End: %f",
                initial_score,
                final_score);
  }
  return -final_score;
}

static double
bfgs_params(model_params_t &                                    initial_params,
            size_t                                              partition_index,
            double                                              p_min,
            double                                              p_max,
            double                                              epsilon,
            double                                              pgtol,
            double                                              factor,
            std::function<double()>                             compute_lh,
            std::function<void(size_t, const model_params_t &)> set_func) {
  int task     = START;
  int n_params = static_cast<int>(initial_params.size());
  set_func(partition_index, initial_params);
  double              score         = compute_lh();
  double              initial_score = score;
  int                 csave;
  std::vector<double> gradient(static_cast<size_t>(n_params), 0.0);
  size_t              max_corrections = 20;
  std::vector<double> wa((2 * max_corrections + 5)
                                 * static_cast<size_t>(n_params)
                             + 12 * max_corrections * (max_corrections + 1),
                         0.0);
  std::vector<int>    iwa(3 * static_cast<size_t>(n_params), 0);

  std::vector<double> parameters(initial_params);
  std::vector<double> param_min(static_cast<size_t>(n_params), p_min);
  std::vector<double> param_max(static_cast<size_t>(n_params), p_max);
  logical             lsave[4];
  int                 isave[44];
  double              dsave[29];
  std::vector<int>    bound_type(static_cast<size_t>(n_params), 2);
  int                 iprint = -1;
  size_t              iters  = 0;

  while (iters < 500) {
    debug_string(EMIT_LEVEL_DEBUG, "Running a bfgs iteration");

    setulb(&n_params,
           (int *)&max_corrections,
           parameters.data(),
           param_min.data(),
           param_max.data(),
           bound_type.data(),
           &score,
           gradient.data(),
           &factor,
           &pgtol,
           wa.data(),
           iwa.data(),
           &task,
           &iprint,
           &csave,
           lsave,
           isave,
           dsave);

    debug_print(EMIT_LEVEL_INFO, "BFGS Iter: %lu Score: %.5f", iters, -score);

    set_func(partition_index, parameters);
    score = compute_lh();
    if (IS_FG(task)) {
      for (size_t i = 0; i < static_cast<size_t>(n_params); ++i) {
        double h = epsilon * fabs(parameters[i]);
        if (h < epsilon) { h = epsilon; }
        double temp = parameters[i];
        parameters[i] += h;
        set_func(partition_index, parameters);
        double dlh = compute_lh();
        assert_string(std::isfinite(dlh), "dlh is not finite");
        gradient[i] = (dlh - score) / h;
        assert_string(std::isfinite(gradient[i]), "gradient is not finite");
        parameters[i] = temp;
      }
    } else if (task != NEW_X) {
      break;
    }
    iters++;
  }
  set_func(partition_index, parameters);
  score = compute_lh();
  // assert_string(initial_score >= score, "Failed to improve the
  // likelihood");
  if (initial_score >= score) {
    std::swap(parameters, initial_params);
  } else {
    debug_print(
        EMIT_LEVEL_WARNING,
        "Failed to improve the likelihood after BFGS. Start: %f, End: %f",
        initial_score,
        score);
  }
  return score;
}

double model_t::bfgs_rates(model_params_t &                    initial_rates,
                           const std::vector<pll_operation_t> &ops,
                           const std::vector<unsigned int>     pmatrix_indices,
                           const std::vector<double>           branch_lengths,
                           size_t                              partition_index,
                           double                              pgtol,
                           double                              factor) {
  constexpr double p_min   = 1e-4;
  constexpr double p_max   = 1e4;
  constexpr double epsilon = 1e-4;

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs params");
  return bfgs_params(
      initial_rates,
      partition_index,
      p_min,
      p_max,
      epsilon,
      pgtol,
      factor,
      [&, this]() -> double {
        return -this->compute_lh_partition(
            partition_index, ops, pmatrix_indices, branch_lengths);
      },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_subst_rates(pi, mp);
      });
}

double model_t::gd_rates(model_params_t &       initial_rates,
                         const root_location_t &rl,
                         size_t                 partition_index) {
  constexpr double p_min   = 1e-4;
  constexpr double p_max   = 1e4;
  constexpr double epsilon = 1e-4;
  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs params");
  return gd_params(
      initial_rates,
      partition_index,
      p_min,
      p_max,
      epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_subst_rates(pi, mp);
      });
}

double model_t::bfgs_freqs(model_params_t &                    initial_freqs,
                           const std::vector<pll_operation_t> &ops,
                           const std::vector<unsigned int>     pmatrix_indices,
                           const std::vector<double>           branch_lengths,
                           size_t                              partition_index,
                           double                              pgtol,
                           double                              factor) {
  constexpr double p_min   = 1e-4;
  constexpr double p_max   = 1.0 - 1e-4 * 3;
  constexpr double epsilon = 1e-4;
  model_params_t   constrained_freqs(initial_freqs.begin(),
                                   initial_freqs.end() - 1);

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs freqs");
  double lh = bfgs_params(
      initial_freqs,
      partition_index,
      p_min,
      p_max,
      epsilon,
      pgtol,
      factor,
      [&, this]() -> double {
        return -this->compute_lh_partition(
            partition_index, ops, pmatrix_indices, branch_lengths);
      },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_freqs_all_free(pi, mp);
      });

  return lh;
}

double model_t::gd_freqs(model_params_t &       initial_freqs,
                         const root_location_t &rl,
                         size_t                 partition_index) {
  constexpr double p_min   = 1e-4;
  constexpr double p_max   = 1.0 - 1e-4 * 3;
  constexpr double epsilon = 1e-4;
  model_params_t   constrained_freqs(initial_freqs.begin(),
                                   initial_freqs.end() - 1);

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs freqs");
  double lh = gd_params(
      initial_freqs,
      partition_index,
      p_min,
      p_max,
      epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_freqs_all_free(pi, mp);
      });

  return lh;
}

double
model_t::bfgs_gamma_rates(model_params_t &                    alpha,
                          const std::vector<pll_operation_t> &ops,
                          const std::vector<unsigned int>     pmatrix_indices,
                          const std::vector<double>           branch_lengths,
                          size_t                              partition_index,
                          double                              pgtol,
                          double                              factor) {
  constexpr double p_min   = 0.2;
  constexpr double p_max   = 10000.0;
  constexpr double epsilon = 1e-4;

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs gamma");
  double lh = bfgs_params(
      alpha,
      partition_index,
      p_min,
      p_max,
      epsilon,
      pgtol,
      factor,
      [&, this]() -> double {
        return -this->compute_lh_partition(
            partition_index, ops, pmatrix_indices, branch_lengths);
      },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_gamma_rates(pi, mp);
      });

  return lh;
}

double model_t::gd_gamma_rates(model_params_t &       alpha,
                               const root_location_t &rl,
                               size_t                 partition_index) {
  constexpr double p_min   = 1e-4;
  constexpr double p_max   = 1.0 - 1e-4 * 3;
  constexpr double epsilon = 1e-4;

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs gamma");
  double lh = gd_params(
      alpha,
      partition_index,
      p_min,
      p_max,
      epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_gamma_rates(pi, mp);
      });

  return lh;
}

double
model_t::bfgs_gamma_weights(model_params_t &                    alpha,
                            const std::vector<pll_operation_t> &ops,
                            const std::vector<unsigned int>     pmatrix_indices,
                            const std::vector<double>           branch_lengths,
                            size_t                              partition_index,
                            double                              pgtol,
                            double                              factor) {
  constexpr double p_min   = 1e-4;
  constexpr double p_max   = 1.0;
  constexpr double epsilon = 1e-4;

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs gamma");
  double lh = bfgs_params(
      alpha,
      partition_index,
      p_min,
      p_max,
      epsilon,
      pgtol,
      factor,
      [&, this]() -> double {
        return -this->compute_lh_partition(
            partition_index, ops, pmatrix_indices, branch_lengths);
      },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_gamma_weights(pi, mp);
      });

  return lh;
}

double model_t::gd_gamma_weights(model_params_t &       alpha,
                                 const root_location_t &rl,
                                 size_t                 partition_index) {
  constexpr double p_min   = 1e-4;
  constexpr double p_max   = 1.0 - 1e-4 * 3;
  constexpr double epsilon = 1e-4;

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs gamma");
  double lh = gd_params(
      alpha,
      partition_index,
      p_min,
      p_max,
      epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_gamma_weights(pi, mp);
      });

  return lh;
}

std::vector<double> model_t::compute_all_root_lh() {
  compute_lh(_tree.roots()[0]);
  std::vector<double> root_lh;
  root_lh.reserve(_tree.roots().size());
  for (auto rl : _tree.roots()) {
    move_root(rl);
    root_lh.push_back(compute_lh_root(rl));
  }
  return root_lh;
}

void model_t::set_subst_rates_uniform() {
  for (size_t i = 0; i < _partitions.size(); ++i) {
    unsigned int   states = _partitions[i]->states;
    unsigned int   params = states * states - states;
    model_params_t mp(params, 1.0 / params);
    set_subst_rates(i, mp);
  }
}

void model_t::set_empirical_freqs() {
  for (size_t i = 0; i < _partitions.size(); ++i) { set_empirical_freqs(i); }
}

void model_t::assign_indicies(const std::vector<size_t> &idx) {
  _assigned_idx = idx;
}

void model_t::assign_indicies(size_t begin, size_t end) {
  assert_string(end >= begin, "We got a negative value for a resize");
  _assigned_idx.resize(end - begin);
  if (_assigned_idx.size() != 0) {
    std::iota(_assigned_idx.begin(), _assigned_idx.end(), begin);
  }
}

void model_t::assign_indicies() {
  _assigned_idx.resize(_tree.root_count());
  std::iota(_assigned_idx.begin(), _assigned_idx.end(), 0);
}

void model_t::assign_indicies(size_t beg, size_t end, std::vector<size_t> idx) {
  _assigned_idx.clear();
  _assigned_idx.reserve(end - beg);
  for (size_t i = beg; i < end; ++i) { _assigned_idx.push_back(idx[i]); }
}

std::pair<size_t, size_t>
model_t::compute_chunk_size_mod(size_t num_tasks) const {
  return compute_chunk_size_mod(_tree.root_count(), num_tasks);
}

std::pair<size_t, size_t>
model_t::compute_chunk_size_mod(size_t root_count, size_t num_tasks) const {
  size_t chunk_size = root_count / num_tasks;
  size_t mod        = root_count % num_tasks;
  return {chunk_size, mod};
}

void model_t::assign_indicies_by_rank_search(size_t        min_roots,
                                             double        root_ratio,
                                             size_t        rank,
                                             size_t        num_tasks,
                                             checkpoint_t &checkpoint) {
  auto   completed_work = checkpoint.completed_indicies();
  auto   shuffled_idx   = shuffle_root_indicies();
  size_t root_count     = std::min(
      std::max(static_cast<size_t>(_tree.root_count() * root_ratio), min_roots),
      _tree.root_count());

  if (root_count < completed_work.size()) {
    throw std::runtime_error{"There are too many results in the checkpoint for "
                             "this search. Is the checkpoint corrupted?"};
  }

  std::sort(completed_work.begin(), completed_work.end());

  size_t work_left = root_count - completed_work.size();
  if (work_left < num_tasks) {
    debug_string(EMIT_LEVEL_WARNING,
                 "The number of searches is less than the number of processes. "
                 "Some nodes will do no work");
  }

  std::vector<size_t> trimmed_idx;
  trimmed_idx.reserve(work_left);

  for (auto i : shuffled_idx) {
    if (!std::binary_search(completed_work.begin(), completed_work.end(), i)) {
      trimmed_idx.push_back(i);
    }
  }

  size_t chunk_size, mod;
  {
    auto res   = compute_chunk_size_mod(work_left, num_tasks);
    chunk_size = res.first;
    mod        = res.second;
  }

  size_t beg = chunk_size * rank + std::min(mod, rank);
  size_t end = chunk_size * (rank + 1) + std::min(mod, (rank + 1));
  assign_indicies(beg, end, trimmed_idx);
}

void model_t::assign_indicies_by_rank_exhaustive(size_t        rank,
                                                 size_t        num_tasks,
                                                 checkpoint_t &checkpoint) {
  auto completed_work = checkpoint.current_progress();
  if (_tree.root_count() < completed_work.size()) {
    throw std::runtime_error{"There are too many results in the checkpoint for "
                             "this tree, are you sure the checkpoint matches?"};
  }
  size_t work_left = _tree.root_count() - completed_work.size();
  debug_print(EMIT_LEVEL_MPI_DEBUG, "work_left: %lu", work_left);
  if (work_left < num_tasks) {
    debug_string(EMIT_LEVEL_WARNING,
                 "The number of searches is less than the number of processes. "
                 "Some nodes will do no work");
  }

  std::sort(completed_work.begin(),
            completed_work.end(),
            [](rd_result_t a, rd_result_t b) { return a.root_id < b.root_id; });

  std::vector<size_t> tmp_idx;
  tmp_idx.reserve(work_left);
  size_t curr_work_idx = 0;
  for (size_t i = 0; i < _tree.root_count(); ++i) {
    if (curr_work_idx < completed_work.size()
        && completed_work[curr_work_idx].root_id == i) {
      ++curr_work_idx;
    } else {
      tmp_idx.push_back(i);
    }
  }

  size_t chunk_size, mod;
  {
    auto res   = compute_chunk_size_mod(work_left, num_tasks);
    chunk_size = res.first;
    mod        = res.second;
  }
  size_t beg = chunk_size * rank + std::min(mod, rank);
  size_t end = chunk_size * (rank + 1) + std::min(mod, (rank + 1));
  assign_indicies(beg, end, tmp_idx);
#ifdef MPI_VERSION
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void model_t::set_model_params(
    const std::vector<partition_parameters_t> &params) {
  for (size_t i = 0; i < params.size(); ++i) {
    set_subst_rates(i, params[i].subst_rates);
    set_freqs(i, params[i].freqs);
    set_gamma_rates(i, params[i].gamma_alpha);
    if (_rate_category_types[i] == rate_category::FREE) {
      set_gamma_weights(i, params[i].gamma_weights);
    }
  }
}

void model_t::optimize_params(std::vector<partition_parameters_t> &params,
                              const root_location_t &              rl,
                              double                               pgtol,
                              double                               factor,
                              bool optimize_gamma) {
  std::vector<pll_operation_t> ops;
  std::vector<unsigned int>    pmatrix_indices;
  std::vector<double>          branch_lengths;

  GENERATE_AND_UNPACK_OPS(_tree, rl, ops, pmatrix_indices, branch_lengths);
#pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < _partitions.size(); ++i) {
    set_subst_rates(i, params[i].subst_rates);
    set_freqs_all_free(i, params[i].freqs);
    set_gamma_rates(i, params[i].gamma_alpha);

    if (_rate_category_types[i] == rate_category::FREE) {
      set_gamma_weights(i, params[i].gamma_weights);
    }

    debug_string(EMIT_LEVEL_INFO, "Optmizing rates");
    bfgs_rates(params[i].subst_rates,
               ops,
               pmatrix_indices,
               branch_lengths,
               i,
               pgtol,
               factor);

    debug_string(EMIT_LEVEL_INFO, "Optimizing freqs");
    bfgs_freqs(params[i].freqs,
               ops,
               pmatrix_indices,
               branch_lengths,
               i,
               pgtol,
               factor);

    if (optimize_gamma) {
      debug_string(EMIT_LEVEL_INFO, "Optimizing gamma rates");
      bfgs_gamma_rates(params[i].gamma_alpha,
                       ops,
                       pmatrix_indices,
                       branch_lengths,
                       i,
                       pgtol,
                       factor);

      if (_rate_category_types[i] == rate_category::FREE) {
        debug_string(EMIT_LEVEL_INFO, "Optimizing gamma weights");
        bfgs_gamma_weights(params[i].gamma_weights,
                           ops,
                           pmatrix_indices,
                           branch_lengths,
                           i,
                           pgtol,
                           factor);
      }
    }
  }
}
