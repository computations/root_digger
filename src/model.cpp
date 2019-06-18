#include "debug.h"
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
                 const msa_t &msa, const model_params_t &freqs) {

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
  if(PLL_STAT(avx2_present)){
    attributes |= PLL_ATTRIB_ARCH_AVX2;
  }
  else if(PLL_STAT(avx_present)){
    attributes |= PLL_ATTRIB_ARCH_AVX;
  }
  else if(PLL_STAT(sse42_present)){
    attributes |= PLL_ATTRIB_ARCH_SSE;
  }

  attributes |= PLL_ATTRIB_NONREV;

  _partition = pll_partition_create(_tree.tip_count(), _tree.branch_count(),
                                    msa.states(), msa.length(), submodels,
                                    _tree.branch_count(), submodels,
                                    _tree.branch_count(), attributes);
  pll_set_subst_params(_partition, 0, rate_parameters.data());

  double rate_cats[submodels] = {0};
  pll_compute_gamma_cats(1, 4, rate_cats, PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates(_partition, rate_cats);
  debug_print("rate: %f", rate_cats[0]);

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

  /* set the frequencies
   *
   * For now we are going to just use emperical frequencies, but in the future,
   * we can do something more clever.
   */
  if (freqs.size() == 0) {
    double *emp_freqs = pllmod_msa_empirical_frequencies(_partition);
    for (size_t i = 0; i < msa.states(); ++i) {
      if (emp_freqs[i] <= 0) {
        free(emp_freqs);
        pll_partition_destroy(_partition);
        throw std::runtime_error("One of the state frequenices is zero while "
                                 "using emperical frequenices");
      }
    }
    pll_set_frequencies(_partition, 0, emp_freqs);
    free(emp_freqs);
  } else {
    for (auto f : freqs) {
      if (f <= 0.0) {
        pll_partition_destroy(_partition);
        throw std::runtime_error("Frequencies with 0 entries are not allowed");
      }
    }
    pll_set_frequencies(_partition, 0, freqs.data());
  }
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

  unsigned int params[] = {0};

  int result =
      pll_update_prob_matrices(_partition, params, pmatrix_indices.data(),
                               branch_lengths.data(), pmatrix_indices.size());

  if (result == PLL_FAILURE) {
    throw std::runtime_error(pll_errmsg);
  }
  pll_update_partials(_partition, &op, 1);
  std::vector<double> persite(_partition->sites);
  double lh = pll_compute_root_loglikelihood(_partition, _tree.root_clv_index(),
                                             _tree.root_scaler_index(), params,
                                             persite.data());
  if (std::isnan(lh)) {
    throw std::runtime_error("lh at root is not a number: " +
                             std::to_string(lh));
  }
  return lh;
}

/*
 * Use a tangent method to compute the derivative
 * TODO: Allow for the return of the likelihood at root, so that we can save an
 * evalutation.
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

std::pair<root_location_t, double> model_t::bisect(const root_location_t &beg,
                                                   dlh_t d_beg,
                                                   const root_location_t &end,
                                                   dlh_t d_end, double atol,
                                                   size_t depth = 0) {
  root_location_t midpoint{beg};
  midpoint.brlen_ratio = (beg.brlen_ratio + end.brlen_ratio) / 2;

  auto d_midpoint = compute_dlh(midpoint);

  if (depth > 64) {
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

/* Find the optimum for the ratio via the bisection method */
root_location_t model_t::optimize_alpha(const root_location_t &root) {
  constexpr double ATOL = 1e-7;
  double lh = compute_lh(root);
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

  if (fabs(d_beg.dlh) < ATOL || fabs(d_end.dlh) < ATOL) {
    debug_string("one of the endpoints is sufficient");
    return best_endpoint;
  }

  if ((d_beg.dlh < 0.0 && d_end.dlh > 0.0) ||
      (d_beg.dlh > 0.0 && d_end.dlh < 0.0)) {
    auto mid = bisect(beg, d_beg, end, d_end, ATOL);
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
      if (fabs(d_midpoint.dlh) < ATOL) {
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
        auto r1 = bisect(beg, d_beg, midpoint_root, d_midpoint, ATOL);
        auto r2 = bisect(midpoint_root, d_midpoint, end, d_end, ATOL);
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

root_location_t model_t::optimize_root_location() {
  std::pair<root_location_t, double> best;
  best.second = -INFINITY;
  for (size_t i = 0; i < _tree.root_count(); ++i) {
    debug_print("working rl: %s", _tree.root_location(i).label().c_str());
    root_location_t rl = optimize_alpha(_tree.root_location(i));
    debug_print("alpha: %f", rl.brlen_ratio);
    double rl_lh = compute_lh(rl);
    debug_print("rl_lh: %f", rl_lh);
    if (rl_lh > best.second) {
      best.first = rl;
      best.second = rl_lh;
    }
  }
  debug_print("finished with lh: %f", best.second);
  return best.first;
}

const rooted_tree_t &model_t::rooted_tree(const root_location_t &root) {
  _tree.root_by(root);
  return _tree;
}
