#include "model.hpp"
#include <cmath>
#include <fstream>
#include <string>
#include <unordered_map>
extern "C" {
#include <libpll/pll_msa.h>
}

#include <iostream>

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
  attributes |= PLL_ATTRIB_ARCH_AVX;
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
  std::cout << "[model_t] rate: " << rate_cats[0] << std::endl;

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

  /*
  std::cout << "[compute_lh] updating partials with these pmatricies"
            << std::endl;
  for (size_t i = 0; i < pmatrix_indices.size(); ++i) {
    std::cout << "brlen: " << branch_lengths[i]
              << ", index: " << pmatrix_indices[i] << std::endl;
    pll_show_pmatrix(_partition, pmatrix_indices[i], 10);
  }
  */

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
    std::cout << "[compute_dlh] both evals are inf, returning 0 for derivative"
              << std::endl;
    return {fx, 0};
  }
  double dlh = (fxh - fx) / EPSILON;
  std::cout << "[compute_dlh] dlh: " << dlh * sign << ", fx: " << fx
            << ", fxh: " << fxh << std::endl;
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
    std::cout << "[bisect] depth too low, returning midpoint with ratio: "
              << midpoint.brlen_ratio << ", lh: " << compute_lh_root(midpoint)
              << std::endl;
    return {midpoint, d_midpoint.lh};
  }

  /* case 1: d_midpoint is within tolerance of 0.0
   * We found a root, and should return
   */

  if (fabs(d_midpoint.dlh) < atol) {
    std::cout << "[bisect] case 1" << std::endl;
    return {midpoint, d_midpoint.lh};
  }

  /* case 2: midpoint is opposite sign of both beg and end
   * There are at least 2 roots, so we recurse on both sides, and return the one
   * with the best LH
   */
  if ((d_beg.dlh > 0.0 && d_end.dlh > 0.0 && d_midpoint.dlh < 0.0) ||
      (d_beg.dlh < 0.0 && d_end.dlh < 0.0 && d_midpoint.dlh > 0.0)) {
    std::cout << "[bisect] case 2" << std::endl;
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
    std::cout << "[bisect] case 3, d_beg: " << d_beg.dlh
              << ", d_midpoint: " << d_midpoint.dlh << ", d_end: " << d_end.dlh
              << ", alpha: " << midpoint.brlen_ratio << std::endl;
    return bisect(beg, d_beg, midpoint, d_midpoint, atol, depth + 1);
  }

  /* case 4: beg and midpoint share a sign, while end has the opposite sign
   * In this case, there is a root between midpoint and end
   */
  if ((d_beg.dlh < 0.0 && d_midpoint.dlh < 0.0 && d_end.dlh > 0.0) ||
      (d_beg.dlh > 0.0 && d_midpoint.dlh > 0.0 && d_end.dlh < 0.0)) {
    std::cout << "[bisect] case 4" << std::endl;
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
    std::cout << "[optimize_alpha] beg end point is best" << std::endl;
  } else {
    std::cout << "[optimize_alpha] end end point is best" << std::endl;
  }

  std::cout << "[optimize_alpha] lh_endpoint.dlh: " << lh_best_endpoint.dlh
            << std::endl;
  std::cout << "[optimize_alpha] d_beg.dlh: " << d_beg.dlh
            << " d_end.dlh: " << d_end.dlh << std::endl;

  if (fabs(d_beg.dlh) < ATOL || fabs(d_end.dlh) < ATOL) {
    std::cout << "[optimize_alpha] one of the endpoints is sufficient"
              << std::endl;
    return best_endpoint;
  }

  if ((d_beg.dlh < 0.0 && d_end.dlh > 0.0) ||
      (d_beg.dlh > 0.0 && d_end.dlh < 0.0)) {
    auto mid = bisect(beg, d_beg, end, d_end, ATOL);
    std::cout << "[optimize_alpha] mid lh: " << mid.second
              << " end lh: " << lh_best_endpoint.lh << std::endl;
    return lh_best_endpoint.lh > mid.second ? best_endpoint : mid.first;
  }

  /* if we have gotten here, then we must have an "even" function so, grid
   * search to find a midpoint which has an opposite sign value and use that to
   * do bisection with
   */

  bool beg_end_pos = d_beg.dlh > 0.0 && d_end.dlh > 0.0;

  for (size_t midpoints = 2; midpoints <= 64; midpoints *= 2) {
    for (size_t midpoint = 1; midpoint <= midpoints; ++midpoint) {
      if (midpoint % 2 == 0)
        continue;

      double alpha = 1.0 / (double)midpoints * midpoint;
      std::cout << "[optimize_alpha] alpha: " << alpha << std::endl;
      root_location_t midpoint_root{beg};
      midpoint_root.brlen_ratio = alpha;
      auto d_midpoint = compute_dlh(midpoint_root);
      std::cout << "[optimize_alpha] d_midpoint.dlh: " << d_midpoint.dlh
                << std::endl;
      if (fabs(d_midpoint.dlh) < ATOL) {
        return midpoint_root;
      }
      if ((beg_end_pos && d_midpoint.dlh < 0.0) ||
          (!beg_end_pos && d_midpoint.dlh > 0.0)) {
        /*
         * we have a midpoint, so now we need to figure out if it is a min or a
         * maximum.
         */
        auto r1 = bisect(beg, d_beg, midpoint_root, d_midpoint, ATOL);
        auto r2 = bisect(midpoint_root, d_midpoint, end, d_end, ATOL);
        std::cout << "[optimize_alpha] r1 lh: " << r1.second
                  << " r2 lh: " << r2.second << std::endl;
        if (r1.second < r2.second) {
          return lh_best_endpoint.lh >= r2.second ? best_endpoint : r2.first;
        }
        return lh_best_endpoint.lh >= r1.second ? best_endpoint : r1.first;
      }
    }
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
    std::cout << "[optimize_root_location] working rl: "
              << _tree.root_location(i).label() << std::endl;
    root_location_t rl = optimize_alpha(_tree.root_location(i));
    std::cout << "[optimize_root_location] alpha: " << rl.brlen_ratio
              << std::endl;
    double rl_lh = compute_lh(rl);
    std::cout << "[optimize_root_location] rl_lh: " << rl_lh << std::endl;
    if (rl_lh > best.second) {
      best.first = rl;
      best.second = rl_lh;
    }
  }
  std::cout << "[optimize_root_location] finished with lh: " << best.second
            << std::endl;
  return best.first;
}

const rooted_tree_t &model_t::rooted_tree(const root_location_t &root) {
  _tree.root_by(root);
  return _tree;
}
