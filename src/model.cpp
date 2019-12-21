#include "model.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
extern "C" {
#include <lbfgsb.h>
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

model_params_t random_params(size_t size, uint64_t seed) {
  model_params_t mp(size);
  std::minstd_rand engine(seed);
  std::uniform_real_distribution<> dist(1e-4, 1.0);
  for (auto &f : mp) {
    f = dist(engine);
  }
  return mp;
}

std::vector<double> sample_dirichlet(std::minstd_rand engine, double alpha,
                                     double beta, size_t size) {
  std::vector<double> rand_sample(size);
  static std::gamma_distribution<double> gd(alpha * (beta / (size * beta)),
                                            1.0);
  double sum = 0.0;
  for (auto &f : rand_sample) {
    f = gd(engine);
    sum += f;
  }
  for (auto &f : rand_sample) {
    f /= sum;
  }
  return rand_sample;
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

void model_t::set_gamma_rates(size_t p_index) {
  pll_compute_gamma_cats(1.0, _rate_weights[p_index].size(),
                         _rate_weights[p_index].data(), PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates(_partitions[p_index], _rate_weights[p_index].data());
}

void model_t::set_gamma_rates(size_t p_index, double alpha) {
  pll_compute_gamma_cats(alpha, _rate_weights[p_index].size(),
                         _rate_weights[p_index].data(), PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates(_partitions[p_index], _rate_weights[p_index].data());
}

void model_t::update_invariant_sites(size_t p_index) {
  pll_update_invariant_sites(_partitions[p_index]);
}

void model_t::set_tip_states(size_t p_index, const msa_t &msa) {
  /* make a label map */
  auto label_map = _tree.label_map();

  /* use the label map to assign tip states in the partition */

  for (int i = 0; i < msa.count(); ++i) {
    try {
      auto result =
          pll_set_tip_states(_partitions[p_index], label_map.at(msa.label(i)),
                             msa.map(), msa.sequence(i));
      if (result == PLL_FAILURE) {
        throw std::runtime_error("failed to set tip " + std::to_string(i));
      }
    } catch (const std::exception &e) {
      throw std::runtime_error(std::string("Could not find taxa ") +
                               msa.label(i) + " in tree");
    }
  }

  /* set pattern weights */
  pll_set_pattern_weights(_partitions[p_index], msa.weights());
}

void model_t::set_empirical_freqs(size_t p_index) {
  pll_partition_t *partition = _partitions[p_index];
  double *emp_freqs = pllmod_msa_empirical_frequencies(partition);
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
  for (auto p : freqs) {
    sum += p;
  }
  for (auto &f : freqs) {
    f /= sum;
  }
  set_freqs(p_index, freqs);
}

model_t::model_t(rooted_tree_t tree, const std::vector<msa_t> &msas,
                 const std::vector<size_t> &rate_cats, uint64_t seed,
                 bool early_stop)
    : _seed{seed}, _early_stop{early_stop} {

  _random_engine = std::minstd_rand(_seed);
  _tree = std::move(tree);
  for (auto rc : rate_cats) {
    _rate_weights.emplace_back(rc, 0.0);
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
    _partitions.push_back(pll_partition_create(
        _tree.tip_count(), _tree.branch_count(), msa.states(), msa.length(),
        _submodels, _tree.branch_count(), _rate_weights[partition_index].size(),
        _tree.branch_count(), attributes));
    _partition_weights.push_back(msa.total_weight());
    set_gamma_rates(partition_index);
    total_weight += msa.total_weight();
  }

  for (auto &&pw : _partition_weights) {
    pw /= (double)total_weight;
  }
}

model_t::~model_t() {
  for (auto p : _partitions) {
    if (p)
      pll_partition_destroy(p);
  }
}

double model_t::compute_lh(const root_location_t &root_location) {
  std::vector<pll_operation_t> ops;
  std::vector<unsigned int> pmatrix_indices;
  std::vector<double> branch_lengths;

  GENERATE_AND_UNPACK_OPS(_tree, root_location, ops, pmatrix_indices,
                          branch_lengths);

  double lh = 0.0;

  for (size_t i = 0; i < _partitions.size(); ++i) {
    auto &partition = _partitions[i];
    int result = pll_update_prob_matrices(
        partition, _param_indicies[i].data(), pmatrix_indices.data(),
        branch_lengths.data(), pmatrix_indices.size());

    if (result == PLL_FAILURE) {
      throw std::runtime_error(pll_errmsg);
    }

    pll_update_partials(partition, ops.data(), ops.size());

    lh += pll_compute_root_loglikelihood(partition, _tree.root_clv_index(),
                                         _tree.root_scaler_index(),
                                         _param_indicies[i].data(), nullptr) *
          _partition_weights[i];
  }
  return lh;
}

double model_t::compute_lh_root(const root_location_t &root) {
  pll_operation_t op;
  std::vector<unsigned int> pmatrix_indices;
  std::vector<double> branch_lengths;

  {
    auto result = _tree.generate_derivative_operations(root);
    op = std::move(std::get<0>(result));
    pmatrix_indices = std::move(std::get<1>(result));
    branch_lengths = std::move(std::get<2>(result));
  }

  unsigned int params[] = {0, 0, 0, 0};
  double lh = 0.0;

  for (size_t i = 0; i < _partitions.size(); ++i) {
    auto &partition = _partitions[i];
    int result =
        pll_update_prob_matrices(partition, params, pmatrix_indices.data(),
                                 branch_lengths.data(), pmatrix_indices.size());

    if (result == PLL_FAILURE) {
      throw std::runtime_error(pll_errmsg);
    }
    pll_update_partials(partition, &op, 1);
    lh += pll_compute_root_loglikelihood(partition, _tree.root_clv_index(),
                                         _tree.root_scaler_index(), params,
                                         nullptr) *
          _partition_weights[i];
  }
  if (std::isnan(lh)) {
    throw std::runtime_error("lh at root is not a number: " +
                             std::to_string(lh));
  }
  return lh;
}

/*
 * Use a secant method to compute the derivative
 */
dlh_t model_t::compute_dlh(const root_location_t &root) {

  dlh_t ret;
  constexpr double EPSILON = 1e-8;
  root_location_t root_prime{root};
  root_prime.brlen_ratio += EPSILON;
  double sign = 1.0;
  if (root_prime.brlen_ratio >= 1.0) {
    root_prime.brlen_ratio = root.brlen_ratio - EPSILON;
    sign = -1.0;
  }

  double fx = compute_lh_root(root);
  ret.lh = fx;

  if (std::isnan(fx)) {
    throw std::runtime_error("fx is not finite when computing derivative: " +
                             std::to_string(root.edge->length));
  }
  double fxh = compute_lh_root(root_prime);
  if (std::isnan(fxh)) {
    throw std::runtime_error("fxh is not finite when computing derivative: " +
                             std::to_string(root_prime.edge->length));
  }

  if (std::isnan(fxh)) {
    throw std::runtime_error("fxh is not finite when computing derivative: " +
                             std::to_string(root_prime.edge->length));
  }
  if (std::isinf(fxh) && std::isinf(fx)) {
    debug_string(EMIT_LEVEL_DEBUG,
                 "Both evals are -inf, returning a 0 derivative");
    return {fx, 0};
  }
  double dlh = (fxh - fx) / EPSILON;
  debug_print(EMIT_LEVEL_DEBUG, "dlh: %f, fx: %f, fxh: %f", dlh * sign, fx,
              fxh);
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
                                                   dlh_t d_beg,
                                                   const root_location_t &end,
                                                   dlh_t d_end, double atol,
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
        midpoint.brlen_ratio, d_midpoint.lh);
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
   * There are at least 2 roots, so we recurse on both sides, and return the one
   * with the best LH
   */
  if ((d_beg.dlh > 0.0 && d_end.dlh > 0.0 && d_midpoint.dlh < 0.0) ||
      (d_beg.dlh < 0.0 && d_end.dlh < 0.0 && d_midpoint.dlh > 0.0)) {
    debug_string(EMIT_LEVEL_DEBUG, "case 2");
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
  if ((d_beg.dlh < 0.0 && d_midpoint.dlh > 0.0 && d_end.dlh > 0.0) ||
      (d_beg.dlh > 0.0 && d_midpoint.dlh < 0.0 && d_end.dlh < 0.0)) {
    debug_print(
        EMIT_LEVEL_DEBUG,
        "case 3, d_beg.dlh: %f, d_midpoint.dlh: %f, d_end:%f, alpha: %f",
        d_beg.lh, d_midpoint.dlh, d_end.dlh, midpoint.brlen_ratio);
    return bisect(beg, d_beg, midpoint, d_midpoint, atol, depth + 1);
  }

  /* case 4: beg and midpoint share a sign, while end has the opposite sign
   * In this case, there is a root between midpoint and end
   */
  if ((d_beg.dlh < 0.0 && d_midpoint.dlh < 0.0 && d_end.dlh > 0.0) ||
      (d_beg.dlh > 0.0 && d_midpoint.dlh > 0.0 && d_end.dlh < 0.0)) {
    debug_string(EMIT_LEVEL_DEBUG, "case 4");
    return bisect(midpoint, d_midpoint, end, d_end, atol, depth + 1);
  }

  /* case 5: something went wrong */
  throw std::runtime_error("Bisection failed to converge with interval : " +
                           std::to_string(d_beg.dlh) + ", " +
                           std::to_string(d_end.dlh) +
                           ", midpoint.dlh: " + std::to_string(d_midpoint.dlh));
}

std::pair<root_location_t, double> model_t::brents(root_location_t beg,
                                                   dlh_t d_beg,
                                                   root_location_t end,
                                                   dlh_t d_end, double atol) {

  assert_string(d_beg.dlh * d_end.dlh < 0,
                "Brents called with endpoints which don't bracket");

  root_location_t midpoint{end};
  auto d_midpoint = d_end;
  double e;
  double d;
  d = e = end.brlen_ratio - beg.brlen_ratio;

  for (size_t i = 0; i < 64; ++i) {
    if (d_end.dlh * d_midpoint.dlh > 0.0) {
      midpoint = beg;
      d_midpoint = d_beg;
      d = e = end.brlen_ratio - beg.brlen_ratio;
    }
    if (abs(d_end.dlh) < abs(d_midpoint.dlh)) {
      beg = end;
      end = midpoint;
      midpoint = beg;
      d_beg = d_end;
      d_end = d_midpoint;
      d_midpoint = d_beg;
    }

    double tol =
        2.0 * abs(end.brlen_ratio) * std::numeric_limits<double>::epsilon() +
        0.5 * atol;
    double e_tol = 0.5 * (midpoint.brlen_ratio - end.brlen_ratio);
    if (abs(e_tol) <= tol || d_end.dlh == 0)
      return {end, d_end.lh};
    if (abs(e) >= tol && abs(d_beg.dlh) > abs(d_end.dlh)) {
      double s = d_end.dlh / d_beg.dlh;
      double p, q;
      if (beg.brlen_ratio == midpoint.brlen_ratio) {
        p = 2.0 * e_tol * s;
        q = 1.0 - s;
      } else {
        q = d_beg.dlh / d_midpoint.dlh;
        double r = d_end.dlh / d_midpoint.dlh;
        p = s * (2.0 * e_tol * q * (q - r) -
                 (end.brlen_ratio - beg.brlen_ratio) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0)
        q = -q;
      p = abs(p);
      double min1 = 3.0 * e_tol * q - abs(e_tol * q);
      double min2 = abs(e * q);
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
    beg = end;
    d_beg = d_end;
    if (abs(d) > tol)
      end.brlen_ratio += d;
    else {
      end.brlen_ratio += e_tol >= 0.0 ? tol : -tol;
    }
    d_end = compute_dlh(end);
  }
  throw std::runtime_error("Brents method failed to converge");
}

/* Find the optimum for the ratio via brents method */
root_location_t model_t::optimize_alpha(const root_location_t &root,
                                        double atol) {
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
        "Initial derivatives failed when optimizing alpha: " +
        std::to_string(root.edge->length));
  }

  root_location_t best_endpoint = d_beg.lh >= d_end.lh ? beg : end;
  auto lh_best_endpoint = d_beg.lh >= d_end.lh ? d_beg : d_end;
  if (d_beg.lh >= d_end.lh) {
    debug_string(EMIT_LEVEL_DEBUG, "beg endpoint is best");
  } else {
    debug_string(EMIT_LEVEL_DEBUG, "end endpoint is best");
  }

  debug_print(EMIT_LEVEL_DEBUG, "lh_endpoint.dlh: %f", lh_best_endpoint.dlh);
  debug_print(EMIT_LEVEL_DEBUG, "d_beg.dlh: %f, d_end.dlh: %f", d_beg.dlh,
              d_end.dlh);

  if (fabs(d_beg.dlh) < atol || fabs(d_end.dlh) < atol) {
    debug_string(EMIT_LEVEL_DEBUG, "one of the endpoints is sufficient");
    return best_endpoint;
  }

  if ((d_beg.dlh < 0.0 && d_end.dlh > 0.0) ||
      (d_beg.dlh > 0.0 && d_end.dlh < 0.0)) {
    auto mid = brents(beg, d_beg, end, d_end, atol);
    debug_print(EMIT_LEVEL_DEBUG, "mid lh: %f, end lh: %f", mid.second,
                lh_best_endpoint.lh);
    return lh_best_endpoint.lh > mid.second ? best_endpoint : mid.first;
  }

  /* if we have gotten here, then we must have an "even" function so, grid
   * search to find a midpoint which has an opposite sign value and use that to
   * do bisection with
   */

  bool beg_end_pos = d_beg.dlh > 0.0 && d_end.dlh > 0.0;
  dlh_t best_midpoint_lh = {-INFINITY, 0};
  root_location_t best_midpoint;
  bool found_midpoint = false;

  for (size_t midpoints = 2; midpoints <= 32; midpoints *= 2) {
    for (size_t midpoint = 1; midpoint <= midpoints; ++midpoint) {
      if (midpoint % 2 == 0)
        continue;

      double alpha = 1.0 / (double)midpoints * midpoint;
      debug_print(EMIT_LEVEL_DEBUG, "alpha: %f", alpha);
      root_location_t midpoint_root{beg};
      midpoint_root.brlen_ratio = alpha;
      auto d_midpoint = compute_dlh(midpoint_root);
      debug_print(EMIT_LEVEL_DEBUG, "d_midpoint.dlh: %f", d_midpoint.dlh);
      if (fabs(d_midpoint.dlh) < atol) {
        if (best_midpoint_lh.lh < d_midpoint.lh) {
          best_midpoint_lh = d_midpoint;
          best_midpoint = midpoint_root;
          found_midpoint = true;
        }
      }
      if ((beg_end_pos && d_midpoint.dlh < 0.0) ||
          (!beg_end_pos && d_midpoint.dlh > 0.0)) {
        /*
         * we have a midpoint, so now we need to figure out if it is a min or a
         * maximum.
         */
        auto r1 = brents(beg, d_beg, midpoint_root, d_midpoint, atol);
        auto r2 = brents(midpoint_root, d_midpoint, end, d_end, atol);
        debug_print(EMIT_LEVEL_DEBUG, "r1 lh: %f, r2 lh: %f", r1.second,
                    r2.second);
        if (lh_best_endpoint.lh < best_midpoint_lh.lh) {
          lh_best_endpoint = best_midpoint_lh;
          best_endpoint = best_midpoint;
        }
        if (r1.second < r2.second) {
          return lh_best_endpoint.lh >= r2.second ? best_endpoint : r2.first;
        }
        return lh_best_endpoint.lh >= r1.second ? best_endpoint : r1.first;
      }
    }
  }

  if (found_midpoint) {
    return best_midpoint;
  }

  /*
   * If we got here, we can just return the best of beg or end, since it seems
   * like there is no root. If both beg and end are positive, then that means
   * the liklihood _increases_ as we increase alpha, so the best liklihood is at
   * the end. Likewise, if both are negative, then the best endpoint is at the
   * begining.
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
  best.second = -INFINITY;
  _tree.root_by(_tree.root_location(0));

  /* start by making a list of "good" roots, with the current model*/

  auto sorted_roots = suggest_roots(min_roots, root_ratio);
  for (auto &sr : sorted_roots) {
    auto &rl = sr.first;
    debug_print(EMIT_LEVEL_DEBUG, "working rl: %s", rl.label().c_str());
    move_root(rl);
    rl = optimize_alpha(rl, 1e-14);
    debug_print(EMIT_LEVEL_DEBUG, "alpha: %f", rl.brlen_ratio);
    double rl_lh = compute_lh_root(rl);
    debug_print(EMIT_LEVEL_DEBUG, "rl_lh: %f", rl_lh);
    if (rl_lh > best.second) {
      best.first = rl;
      best.second = rl_lh;
    }
  }
  debug_print(EMIT_LEVEL_DEBUG, "finished with lh: %f", best.second);
  return best;
}

void model_t::move_root(const root_location_t &root_location) {
  std::vector<pll_operation_t> ops;
  std::vector<unsigned int> pmatrix_indices;
  std::vector<double> branch_lengths;

  GENERATE_AND_UNPACK_OPS(_tree, root_location, ops, pmatrix_indices,
                          branch_lengths);

  for (size_t i = 0; i < _partitions.size(); ++i) {
    auto &partition = _partitions[i];
    int result = pll_update_prob_matrices(
        partition, _param_indicies[i].data(), pmatrix_indices.data(),
        branch_lengths.data(), pmatrix_indices.size());

    if (result == PLL_FAILURE) {
      throw std::runtime_error(pll_errmsg);
    }

    pll_update_partials(partition, ops.data(), ops.size());
  }
}

std::vector<std::pair<root_location_t, double>> model_t::suggest_roots() {
  return suggest_roots(1, 0.0);
}

std::vector<std::pair<root_location_t, double>>
model_t::suggest_roots(size_t min, double ratio) {
  std::vector<std::pair<root_location_t, double>> rl_lhs;
  rl_lhs.reserve(_tree.root_count());
  for (auto rl : _tree.roots()) {
    move_root(rl);
    rl_lhs.push_back(std::make_pair(rl, compute_lh_root(rl)));
  }
  size_t final_size = std::max(static_cast<size_t>(rl_lhs.size() * ratio), min);
  std::partial_sort(
      rl_lhs.begin(), rl_lhs.begin() + final_size, rl_lhs.end(),
      [](std::pair<root_location_t, double> a,
         std::pair<root_location_t, double> b) { return a.second > b.second; });
  rl_lhs.resize(final_size);
  return rl_lhs;
}

std::vector<root_location_t> model_t::suggest_roots_random(size_t min,
                                                           double ratio) {
  std::vector<root_location_t> random_roots;
  size_t root_count =
      std::max(static_cast<size_t>(_tree.root_count() * ratio), min);
  random_roots.reserve(root_count);
  std::uniform_int_distribution<size_t> dis(0, _tree.root_count() - 1);
  for (size_t i = 0; i < root_count; ++i) {
    size_t root_choice_index = dis(_random_engine);
    random_roots.emplace_back(_tree.root_location(root_choice_index));
  }
  return random_roots;
}

/* Optimize the substitution parameters and root location.*/
std::pair<root_location_t, double>
model_t::optimize_all(size_t min_roots, double root_ratio, double atol,
                      double pgtol, double brtol, double factor) {
  double best_lh = -INFINITY;
  root_location_t best_rl;
  std::vector<model_params_t> initial_subst;
  initial_subst.reserve(_partitions.size());

  for (auto p : _partitions) {
    size_t params_size = p->states * p->states - p->states;
    initial_subst.emplace_back(random_params(params_size, _random_engine()));
  }

  std::vector<model_params_t> initial_freqs;
  initial_freqs.reserve(_partitions.size());
  std::vector<double> gamma_alpha;

  for (size_t i = 0; i < _partitions.size(); ++i) {
    initial_freqs.emplace_back(_partitions[i]->frequencies[0],
                               _partitions[i]->frequencies[0] +
                                   _partitions[i]->states);
    gamma_alpha.push_back(1.0);
  }

  set_subst_rates_uniform();
  set_empirical_freqs();
  auto roots = suggest_roots_random(min_roots, root_ratio);
  size_t root_count = roots.size();
  size_t root_index = 0;

  for (auto rl : roots) {
    set_subst_rates_uniform();
    set_empirical_freqs();
    ++root_index;
    debug_print(EMIT_LEVEL_PROGRESS, "Root %lu/%lu", root_index, root_count);

    move_root(rl);
    std::vector<model_params_t> subst_rates;
    std::vector<model_params_t> freqs;

    for (size_t p = 0; p < _partitions.size(); ++p) {
      size_t subst_size = _partitions[p]->states;
      subst_size = subst_size * subst_size - subst_size;
      freqs.push_back(
          model_params_t(_partitions[p]->states, 1.0 / _partitions[p]->states));
      subst_rates.push_back(model_params_t(subst_size, 1.0 / subst_size));
    }

    auto cur_best_rl = rl;
    double cur_best_lh = -INFINITY;

    for (size_t iter = 0; iter < 1e3; ++iter) {
      for (size_t i = 0; i < _partitions.size(); ++i) {
        set_subst_rates(i, subst_rates[i]);
        set_freqs_all_free(i, freqs[i]);
        set_gamma_rates(i, gamma_alpha[i]);
        debug_string(EMIT_LEVEL_INFO, "Optimizing Gamma");
        bfgs_gamma(gamma_alpha[i], rl, i, pgtol, factor);
        debug_string(EMIT_LEVEL_INFO, "Optmizing Rates");
        bfgs_rates(subst_rates[i], rl, i, pgtol, factor);
        debug_string(EMIT_LEVEL_INFO, "Optimizing Freqs");
        bfgs_freqs(freqs[i], rl, i, pgtol, factor);
      }

      if (fabs(compute_lh(rl) - cur_best_lh) < atol) {
        break;
      }

      debug_string(EMIT_LEVEL_INFO, "Optimizing Root Location");
      auto cur = optimize_root_location(min_roots, root_ratio);

      debug_print(EMIT_LEVEL_INFO, "Iteration %lu LH: %.5f", iter, cur.second);

      if (_early_stop) {
        if (rl.edge == cur.first.edge &&
            fabs(rl.brlen_ratio - cur.first.brlen_ratio) < atol) {
          cur_best_rl = cur.first;
          cur_best_lh = cur.second;
          break;
        }
      }
      if (fabs(cur.second - cur_best_lh) < brtol) {
        cur_best_rl = cur.first;
        cur_best_lh = cur.second;
        break;
      }

      if (cur.second > cur_best_lh) {
        cur_best_rl = cur.first;
        cur_best_lh = cur.second;
      }

      rl = cur_best_rl;
    }

    if (cur_best_lh > best_lh) {
      best_rl = cur_best_rl;
      best_lh = cur_best_lh;
    }
  }
  return {best_rl, best_lh};
}

const rooted_tree_t &model_t::rooted_tree(const root_location_t &root) {
  _tree.root_by(root);
  return _tree;
}

const rooted_tree_t &model_t::virtual_rooted_tree(const root_location_t &root) {
  _tree.root_by(root);
  _tree.unroot();
  return _tree;
}

const rooted_tree_t &model_t::unrooted_tree() {
  _tree.unroot();
  return _tree;
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
    size_t params_size = p->states * p->states - p->states;
    std::vector<double> params{p->subst_params[0],
                               p->subst_params[0] + params_size};
    for (size_t i = 0; i < params.size(); ++i) {
      oss << std::to_string(params[i]);
      if (i != params.size() - 1)
        oss << ",";
    }
    oss << "}";
    if (p_index != _partitions.size() - 1)
      oss << ",";
  }
  oss << "}";
  return oss.str();
}

double gd_params(model_params_t &initial_params, size_t partition_index,
                 double p_min, double p_max, double epsilon,
                 std::function<double()> compute_lh,
                 std::function<void(size_t, const model_params_t &)> set_func) {

  size_t iters = 0;
  double atol = 1e-4;
  size_t n_params = initial_params.size();
  std::vector<double> parameters(initial_params);
  std::vector<double> gradient(n_params, 0.0);
  set_func(partition_index, initial_params);
  double initial_score = compute_lh();
  double last_score = INFINITY;
  while (true) {
    set_func(partition_index, parameters);
    double score = compute_lh();
    debug_print(EMIT_LEVEL_INFO, "GD Iter: %lu Score: %.5f", iters, -score);
    if (fabs(last_score - score) < atol) {
      break;
    }
    for (size_t i = 0; i < n_params; ++i) {
      double h = epsilon * fabs(parameters[i]);
      if (h < epsilon) {
        h = epsilon;
      }
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
      auto tmp_params{parameters};
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
  assert_string(initial_score >= final_score, "Failed to make the lh better");
  std::swap(parameters, initial_params);
  return -final_score;
}

double
bfgs_params(model_params_t &initial_params, size_t partition_index,
            double p_min, double p_max, double epsilon, double pgtol,
            double factor, std::function<double()> compute_lh,
            std::function<void(size_t, const model_params_t &)> set_func) {
  int task = START;
  int n_params = static_cast<int>(initial_params.size());
  set_func(partition_index, initial_params);
  double score = compute_lh();
  double initial_score = score;
  int csave;
  std::vector<double> gradient(n_params, 0.0);
  size_t max_corrections = 20;
  std::vector<double> wa((2 * max_corrections + 5) * n_params +
                             12 * max_corrections * (max_corrections + 1),
                         0.0);
  std::vector<int> iwa(3 * n_params, 0);

  std::vector<double> parameters(initial_params);
  std::vector<double> param_min(n_params, p_min);
  std::vector<double> param_max(n_params, p_max);
  logical lsave[4];
  int isave[44];
  double dsave[29];
  std::vector<int> bound_type(n_params, 2);
  int iprint = -1;
  size_t iters = 0;

  while (iters < 500) {
    debug_string(EMIT_LEVEL_DEBUG, "Running a bfgs iteration");

    setulb(&n_params, (int *)&max_corrections, parameters.data(),
           param_min.data(), param_max.data(), bound_type.data(), &score,
           gradient.data(), &factor, &pgtol, wa.data(), iwa.data(), &task,
           &iprint, &csave, lsave, isave, dsave);

    debug_print(EMIT_LEVEL_DEBUG, "BFGS Iter: %lu Score: %.5f", iters, -score);

    set_func(partition_index, parameters);
    score = compute_lh();
    if (IS_FG(task)) {
      for (int i = 0; i < n_params; ++i) {
        double h = epsilon * fabs(parameters[i]);
        if (h < epsilon) {
          h = epsilon;
        }
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
  assert_string(initial_score >= score, "Failed to improve the likelihood");
  if (initial_score >= score) {
    std::swap(parameters, initial_params);
  }
  return score;
}

double model_t::bfgs_rates(model_params_t &initial_rates,
                           const root_location_t &rl, size_t partition_index,
                           double pgtol, double factor) {

  constexpr double p_min = 1e-4;
  constexpr double p_max = 1e4;
  constexpr double epsilon = 1e-4;
  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs params");
  return bfgs_params(
      initial_rates, partition_index, p_min, p_max, epsilon, pgtol, factor,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_subst_rates(pi, mp);
      });
}

double model_t::gd_rates(model_params_t &initial_rates,
                         const root_location_t &rl, size_t partition_index) {

  constexpr double p_min = 1e-4;
  constexpr double p_max = 1e4;
  constexpr double epsilon = 1e-4;
  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs params");
  return gd_params(
      initial_rates, partition_index, p_min, p_max, epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_subst_rates(pi, mp);
      });
}

double model_t::bfgs_freqs(model_params_t &initial_freqs,
                           const root_location_t &rl, size_t partition_index,
                           double pgtol, double factor) {
  constexpr double p_min = 1e-4;
  constexpr double p_max = 1.0 - 1e-4 * 3;
  constexpr double epsilon = 1e-4;
  model_params_t constrained_freqs(initial_freqs.begin(),
                                   initial_freqs.end() - 1);

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs freqs");
  double lh = bfgs_params(
      initial_freqs, partition_index, p_min, p_max, epsilon, pgtol, factor,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_freqs_all_free(pi, mp);
      });

  return lh;
}

double model_t::gd_freqs(model_params_t &initial_freqs,
                         const root_location_t &rl, size_t partition_index) {
  constexpr double p_min = 1e-4;
  constexpr double p_max = 1.0 - 1e-4 * 3;
  constexpr double epsilon = 1e-4;
  model_params_t constrained_freqs(initial_freqs.begin(),
                                   initial_freqs.end() - 1);

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs freqs");
  double lh = gd_params(
      initial_freqs, partition_index, p_min, p_max, epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_freqs_all_free(pi, mp);
      });

  return lh;
}

double model_t::bfgs_gamma(double &initial_alpha, const root_location_t &rl,
                           size_t partition_index, double pgtol,
                           double factor) {
  constexpr double p_min = 0.2;
  constexpr double p_max = 10000.0;
  constexpr double epsilon = 1e-4;
  model_params_t alpha(1);
  alpha[0] = initial_alpha;

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs gamma");
  double lh = bfgs_params(
      alpha, partition_index, p_min, p_max, epsilon, pgtol, factor,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_gamma_rates(pi, mp[0]);
      });
  initial_alpha = alpha[0];

  return lh;
}

double model_t::gd_gamma(double &initial_alpha, const root_location_t &rl,
                         size_t partition_index) {
  constexpr double p_min = 1e-4;
  constexpr double p_max = 1.0 - 1e-4 * 3;
  constexpr double epsilon = 1e-4;
  model_params_t alpha(1);
  alpha[0] = initial_alpha;

  debug_string(EMIT_LEVEL_DEBUG, "doing bfgs gamma");
  double lh = gd_params(
      alpha, partition_index, p_min, p_max, epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_gamma_rates(pi, mp[0]);
      });
  initial_alpha = alpha[0];

  return lh;
}

std::vector<double> model_t::compute_all_root_lh() {
  _tree.root_by(_tree.roots()[0]);
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
    unsigned int states = _partitions[i]->states;
    unsigned int params = states * states - states;
    model_params_t mp(params, 1.0 / params);
    set_subst_rates(i, mp);
  }
}

void model_t::set_empirical_freqs() {
  for (size_t i = 0; i < _partitions.size(); ++i) {
    set_empirical_freqs(i);
  }
}

std::pair<root_location_t, double> model_t::exhaustive_search(double atol,
                                                              double pgtol,
                                                              double brtol,
                                                              double factor) {
  size_t root_index = 0;
  root_location_t best_rl;
  double best_lh = -INFINITY;
  size_t root_count = _tree.root_count();
  std::vector<std::pair<root_location_t, double>> mapped_likelihoods;
  debug_string(EMIT_LEVEL_PROGRESS, "Starting exhaustive search");

  for (auto rl : _tree.roots()) {
    set_subst_rates_uniform();
    set_empirical_freqs();
    ++root_index;

    move_root(rl);
    std::vector<model_params_t> subst_rates;
    std::vector<model_params_t> freqs;
    std::vector<double> gamma_alphas;

    for (size_t p = 0; p < _partitions.size(); ++p) {
      size_t subst_size = _partitions[p]->states;
      subst_size = subst_size * subst_size - subst_size;
      freqs.push_back(
          model_params_t(_partitions[p]->states, 1.0 / _partitions[p]->states));
      subst_rates.push_back(model_params_t(subst_size, 1.0 / subst_size));
      gamma_alphas.push_back(1.0);
    }

    root_location_t cur_best_rl;
    double cur_best_lh = -INFINITY;

    for (size_t iter = 0; iter < 1e3; ++iter) {
      for (size_t i = 0; i < _partitions.size(); ++i) {
        set_subst_rates(i, subst_rates[i]);
        set_freqs_all_free(i, freqs[i]);
        set_gamma_rates(i, gamma_alphas[i]);
        debug_string(EMIT_LEVEL_DEBUG, "Optimizing gamma rates");
        bfgs_gamma(gamma_alphas[i], rl, i, pgtol, factor);
        debug_string(EMIT_LEVEL_DEBUG, "Optmizing rates");
        bfgs_rates(subst_rates[i], rl, i, pgtol, factor);
        debug_string(EMIT_LEVEL_DEBUG, "Optimizing freqs");
        bfgs_freqs(freqs[i], rl, i, pgtol, factor);
      }

      if (fabs(compute_lh(rl) - cur_best_lh) < atol) {
        break;
      }

      debug_string(EMIT_LEVEL_INFO, "Optimizing Root Location");
      auto cur_rl = optimize_alpha(rl, brtol);
      double cur_lh = compute_lh_root(cur_rl);

      debug_print(EMIT_LEVEL_INFO, "Iteration %lu LH: %.5f", iter, cur_lh);
      debug_print(EMIT_LEVEL_INFO, "difference in lh: %.5f",
                  (cur_lh - cur_best_lh));

      if (_early_stop) {
        if (fabs(rl.brlen_ratio - cur_rl.brlen_ratio) < brtol) {
          debug_print(EMIT_LEVEL_DEBUG,
                      "Current BRlen ratio tolerances: %.7f, brtol: %.7f",
                      fabs(rl.brlen_ratio - cur_rl.brlen_ratio), brtol);
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

    debug_print(EMIT_LEVEL_PROGRESS, "Root %lu / %lu", root_index, root_count);

    mapped_likelihoods.emplace_back(cur_best_rl, cur_best_lh);
    if (cur_best_lh > best_lh) {
      best_rl = cur_best_rl;
      best_lh = cur_best_lh;
    }
  }

  double max_lh = -INFINITY;

  for (auto kv : mapped_likelihoods) {
    max_lh = std::max(kv.second, max_lh);
  }

  double total_lh = 0;
  for (auto kv : mapped_likelihoods) {
    total_lh += exp(kv.second - max_lh);
  }
  debug_print(EMIT_LEVEL_DEBUG, "LWR denom: %f, %e", total_lh, total_lh - 1);

  for (auto kv : mapped_likelihoods) {
    double lwr = exp((kv.second - max_lh)) / total_lh;
    _tree.annotate_branch(kv.first, "LWR", std::to_string(lwr));
    _tree.annotate_lh(kv.first, kv.second);
    _tree.annotate_ratio(kv.first, kv.first.brlen_ratio);
  }

  return {best_rl, best_lh};
}
