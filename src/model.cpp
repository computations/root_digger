#include "debug.h"
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
  set_subst_rates(
      p_index,
      random_params(msa.states() * msa.states() - msa.states(), _seed));
}

void model_t::set_subst_rates_random(size_t p_index, size_t states,
                                     int seed_offset) {
  set_subst_rates(p_index,
                  random_params(states * states - states, _seed + seed_offset));
}

void model_t::set_gamma_rates(size_t p_index) {
  constexpr size_t n_rate_cat = 1;
  double rate_cats[n_rate_cat] = {0};
  pll_compute_gamma_cats(1.0, n_rate_cat, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates(_partitions[p_index], rate_cats);
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
      pll_set_tip_states(_partitions[p_index], label_map.at(msa.label(i)),
                         msa.map(), msa.sequence(i));
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
  double sum = 1.0;
  for (auto p : freqs) {
    sum += p;
  }
  double d = 1.0 / sum;
  for (auto &f : freqs) {
    f *= d;
  }
  freqs.push_back(d);
  set_freqs(p_index, freqs);
}

model_t::model_t(rooted_tree_t tree, const std::vector<msa_t> &msas,
                 uint64_t seed) {

  _seed = seed;
  _random_engine = std::minstd_rand(_seed);
  _tree = std::move(tree);
  _temp_ratio = 0.8;
  _root_opt_frequency = 2.0;

  for (auto &msa : msas) {
    if (!msa.constiency_check(_tree.label_set())) {
      throw std::invalid_argument(
          "Taxa on the tree and in the MSA are inconsistient");
    }
  }

  /*
   * Only one submodel will be used for the time being. If there is desire for
   * more, we can add support for more models..
   */
  constexpr unsigned int submodels = 1;
  constexpr unsigned int n_rate_cat = 1;

  unsigned int attributes = 0;
  if (PLL_STAT(avx2_present)) {
    attributes |= PLL_ATTRIB_ARCH_AVX2;
  } else if (PLL_STAT(avx_present)) {
    attributes |= PLL_ATTRIB_ARCH_AVX;
  } else if (PLL_STAT(sse42_present)) {
    attributes |= PLL_ATTRIB_ARCH_SSE;
  }

  attributes |= PLL_ATTRIB_NONREV;
  // attributes |= PLL_ATTRIB_PATTERN_TIP;

  for (size_t partition_index = 0; partition_index < msas.size();
       ++partition_index) {
    auto &msa = msas[partition_index];
    _partitions.push_back(pll_partition_create(
        _tree.tip_count(), _tree.branch_count(), msa.states(), msa.length(),
        submodels, _tree.branch_count(), n_rate_cat, _tree.branch_count(),
        attributes));
    _partition_weights.push_back(msa.total_weight());
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

  unsigned int params[4] = {0, 0, 0, 0};
  double lh = 0.0;

  for (size_t i = 0; i < _partitions.size(); ++i) {
    auto &partition = _partitions[i];
    int result =
        pll_update_prob_matrices(partition, params, pmatrix_indices.data(),
                                 branch_lengths.data(), pmatrix_indices.size());

    if (result == PLL_FAILURE) {
      throw std::runtime_error(pll_errmsg);
    }

    pll_update_partials(partition, ops.data(), ops.size());

    lh += pll_compute_root_loglikelihood(partition, _tree.root_clv_index(),
                                         _tree.root_scaler_index(), params,
                                         nullptr) *
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
    debug_string("Both evals are -inf, returning a 0 derivative");
    return {fx, 0};
  }
  double dlh = (fxh - fx) / EPSILON;
  debug_print("dlh: %f, fx: %f, fxh: %f", dlh * sign, fx, fxh);
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
        "depth exceeded limit, returning midpoint with ratio: %f, lh: %f",
        midpoint.brlen_ratio, d_midpoint.lh);
    return {midpoint, d_midpoint.lh};
  }

  /* case 1: d_midpoint is within tolerance of 0.0
   * We found a root, and should return
   */

  if (fabs(d_midpoint.dlh) < atol) {
    debug_string("case 1");
    return {midpoint, d_midpoint.lh};
  }

  /* case 2: midpoint is opposite sign of both beg and end
   * There are at least 2 roots, so we recurse on both sides, and return the one
   * with the best LH
   */
  if ((d_beg.dlh > 0.0 && d_end.dlh > 0.0 && d_midpoint.dlh < 0.0) ||
      (d_beg.dlh < 0.0 && d_end.dlh < 0.0 && d_midpoint.dlh > 0.0)) {
    debug_string("case 2");
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
        "case 3, d_beg.dlh: %f, d_midpoint.dlh: %f, d_end:%f, alpha: %f",
        d_beg.lh, d_midpoint.dlh, d_end.dlh, midpoint.brlen_ratio);
    return bisect(beg, d_beg, midpoint, d_midpoint, atol, depth + 1);
  }

  /* case 4: beg and midpoint share a sign, while end has the opposite sign
   * In this case, there is a root between midpoint and end
   */
  if ((d_beg.dlh < 0.0 && d_midpoint.dlh < 0.0 && d_end.dlh > 0.0) ||
      (d_beg.dlh > 0.0 && d_midpoint.dlh > 0.0 && d_end.dlh < 0.0)) {
    debug_string("case 4");
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
    debug_string("beg endpoint is best");
  } else {
    debug_string("end endpoint is best");
  }

  debug_print("lh_endpoint.dlh: %f", lh_best_endpoint.dlh);
  debug_print("d_beg.dlh: %f, d_end.dlh: %f", d_beg.dlh, d_end.dlh);

  if (fabs(d_beg.dlh) < atol || fabs(d_end.dlh) < atol) {
    debug_string("one of the endpoints is sufficient");
    return best_endpoint;
  }

  if ((d_beg.dlh < 0.0 && d_end.dlh > 0.0) ||
      (d_beg.dlh > 0.0 && d_end.dlh < 0.0)) {
    auto mid = brents(beg, d_beg, end, d_end, atol);
    debug_print("mid lh: %f, end lh: %f", mid.second, lh_best_endpoint.lh);
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
      debug_print("alpha: %f", alpha);
      root_location_t midpoint_root{beg};
      midpoint_root.brlen_ratio = alpha;
      auto d_midpoint = compute_dlh(midpoint_root);
      debug_print("d_midpoint.dlh: %f", d_midpoint.dlh);
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
        debug_print("r1 lh: %f, r2 lh: %f", r1.second, r2.second);
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

std::pair<root_location_t, double> model_t::optimize_root_location() {
  std::pair<root_location_t, double> best;
  best.second = -INFINITY;
  _tree.root_by(_tree.root_location(0));

  /* start by making a list of "good" roots, with the current model*/

  std::vector<std::pair<double, root_location_t>> sorted_roots;
  sorted_roots.reserve(_tree.root_count());
  for (const auto &rl : _tree.roots()) {
    _tree.root_by(rl);
    sorted_roots.push_back(std::make_pair(compute_lh(rl), rl));
  }

  std::sort(sorted_roots.begin(), sorted_roots.end(),
            [](const std::pair<double, root_location_t> &a,
               const std::pair<double, root_location_t> &b) {
              return a.first > b.first;
            });

  size_t reduced_size = sorted_roots.size() / 1;
  sorted_roots.resize(std::max(reduced_size, 1ul));
  for (auto &sr : sorted_roots) {
    auto &rl = sr.second;
    debug_print("working rl: %s", rl.label().c_str());
    move_root(rl);
    rl = optimize_alpha(rl, 1e-14);
    debug_print("alpha: %f", rl.brlen_ratio);
    double rl_lh = compute_lh(rl);
    debug_print("rl_lh: %f", rl_lh);
    if (rl_lh > best.second) {
      best.first = rl;
      best.second = rl_lh;
    }
  }
  debug_print("finished with lh: %f", best.second);
  return best;
}

void model_t::move_root(const root_location_t &new_root) {
  compute_lh(new_root);
}

std::vector<std::pair<root_location_t, double>> model_t::suggest_roots() {
  set_subst_rates_uniform();
  set_empirical_freqs();
  std::vector<std::pair<root_location_t, double>> rl_lhs;
  rl_lhs.reserve(_tree.root_count());
  for (auto rl : _tree.roots()) {
    rl_lhs.push_back(std::make_pair(rl, compute_lh(rl)));
  }
  std::sort(
      rl_lhs.begin(), rl_lhs.end(),
      [](std::pair<root_location_t, double> a,
         std::pair<root_location_t, double> b) { return a.second < b.second; });
  rl_lhs.resize(std::min(rl_lhs.size(), 10ul));
  return rl_lhs;
}

/* Optimize the substitution parameters by simulated annealing */
root_location_t model_t::optimize_all() {
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
  for (size_t i = 0; i < _partitions.size(); ++i) {
    initial_freqs.emplace_back(_partitions[i]->frequencies[0],
                               _partitions[i]->frequencies[0] +
                                   _partitions[i]->states);
  }
  size_t iters = 0;
  for (auto rl : _tree.roots()) {
    _tree.root_by(rl);
    for (size_t i = 0; i < _partitions.size(); ++i) {
      size_t subst_size = _partitions[i]->states;
      subst_size = subst_size * subst_size - subst_size;
      auto freqs =
          model_params_t(_partitions[i]->states, 1.0 / _partitions[i]->states);
      auto subst_rates = model_params_t(subst_size, 1.0 / subst_size);
      bfgs_rates(subst_rates, rl, i);
      bfgs_freqs(freqs, rl, i);
    }
    auto cur = optimize_root_location();
    if (cur.second > best_lh) {
      best_rl = cur.first;
      best_lh = cur.second;
    }
    if (iters > 10) {
      break;
    }
    ++iters;
  }
  return best_rl;
}

const rooted_tree_t &model_t::rooted_tree(const root_location_t &root) {
  _tree.root_by(root);
  return _tree;
}

void model_t::initialize_partitions(const std::vector<msa_t> &msa) {
  for (size_t partition_index = 0; partition_index < _partitions.size();
       ++partition_index) {
    set_tip_states(partition_index, msa[partition_index]);
    update_invariant_sites(partition_index);
    set_empirical_freqs(partition_index);
    set_subst_rates_random(partition_index, msa[partition_index]);
    set_gamma_rates(partition_index);
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

void model_t::set_temp_ratio(double t) { _temp_ratio = t; }

void model_t::set_root_opt_frequency(double r) { _root_opt_frequency = r; }

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

double
bfgs_params(model_params_t &initial_params, size_t partition_index,
            double p_min, double p_max, double epsilon,
            std::function<double()> compute_lh,
            std::function<void(size_t, const model_params_t &)> set_func) {
  int task = START;
  size_t n_params = initial_params.size();
  double score = compute_lh();
  double factor = 1e4;
  double pgtol = 1e-4;
  int csave[60];
  std::vector<double> gradient(n_params);
  size_t max_corrections = 5;
  std::vector<double> wa((2 * max_corrections + 5) * n_params +
                         12 * max_corrections * (max_corrections + 1));
  std::vector<int> iwa(3 * n_params);

  std::vector<double> parameters(initial_params);
  std::vector<double> param_min(n_params, p_min);
  std::vector<double> param_max(n_params, p_max);
  logical lsave[4];
  int isave[44];
  double dsave[29];
  std::vector<int> bound_type(n_params, 2);
  int iprint = -1;
  set_func(partition_index, initial_params);
  size_t iters = 0;

  while (iters < 5) {
    debug_string("Running a bfgs iteration");
    setulb((int *)&n_params, (int *)&max_corrections, parameters.data(),
           param_min.data(), param_max.data(), bound_type.data(), &score,
           gradient.data(), &factor, &pgtol, wa.data(), iwa.data(), &task,
           &iprint, csave, lsave, isave, dsave);

    if (IS_FG(task)) {
      score = compute_lh();
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
    } else if (task != NEW_X) {
      break;
    }
    iters++;
  }
  set_func(partition_index, parameters);
  score = compute_lh();
  std::swap(parameters, initial_params);
  return score;
}

double model_t::bfgs_rates(model_params_t &initial_rates,
                           const root_location_t &rl, size_t partition_index) {

  constexpr double p_min = 1e-4;
  constexpr double p_max = 1e4;
  constexpr double epsilon = 1e-7;
  debug_string("doing bfgs params");
  return bfgs_params(
      initial_rates, partition_index, p_min, p_max, epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_subst_rates(pi, mp);
      });
}

double model_t::bfgs_freqs(model_params_t &initial_freqs,
                           const root_location_t &rl, size_t partition_index) {
  constexpr double p_min = 1e-4;
  constexpr double p_max = 1.0 - 1e-4 * 3;
  constexpr double epsilon = 1e-7;
  model_params_t constrained_freqs(initial_freqs.begin(),
                                   initial_freqs.end() - 1);

  double d = *(initial_freqs.end() - 1);
  for (auto &f : constrained_freqs) {
    f /= d;
  }
  debug_string("doing bfgs freqs");
  double lh = bfgs_params(
      constrained_freqs, partition_index, p_min, p_max, epsilon,
      [&, this]() -> double { return -this->compute_lh(rl); },
      [&, this](size_t pi, const model_params_t &mp) -> void {
        this->set_freqs_all_free(pi, mp);
      });

  double sum = 0.0;
  for (size_t i = 0; i < constrained_freqs.size(); ++i) {
    sum += constrained_freqs[i];
  }
  d = 1.0 / sum;
  for (size_t i = 0; i < constrained_freqs.size(); ++i) {
    initial_freqs[i] = constrained_freqs[i] * d;
  }
  *(initial_freqs.end() - 1) = d;
  set_freqs(partition_index, initial_freqs);
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
