#ifndef RD_MODEL_HPP_
#define RD_MODEL_HPP_

extern "C" {
#include <libpll/pll.h>
}
#include "debug.h"
#include "msa.hpp"
#include "tree.hpp"
#include <functional>
#ifdef MPI_BUILD
#include <mpi.h>
#endif
#include <random>
#include <string>
#include <utility>
#include <vector>

typedef std::vector<double> model_params_t;

struct dlh_t {
  double lh;
  double dlh;
};

struct invalid_empirical_frequencies_exception : public std::runtime_error {
  invalid_empirical_frequencies_exception(const char *m)
      : std::runtime_error(m){};
};

std::string read_file_contents(std::ifstream &infile);

double parse_param(std::string::const_iterator begin,
                   std::string::const_iterator end);

model_params_t parse_model_params(const std::string &model_string);

model_params_t parse_model_file(const std::string &model_filename);

model_params_t random_params(size_t size, uint64_t seed);

#ifdef MPI_VERSION
struct exhaustive_mode_results_t {
  int root_id;
  double lh;
  double alpha;
};
constexpr const int exhaustive_mode_results_t_nitems = 3;
#endif

class model_t {
public:
  model_t(
      rooted_tree_t t, const std::vector<msa_t> &msa,
      const std::vector<size_t> &rate_cats,
      const std::vector<rate_category::rate_category_e> &rate_category_types,
      bool invariant_sites, uint64_t seed, bool early_stop);
  model_t(rooted_tree_t t, const std::vector<msa_t> &msa,
          const std::vector<size_t> &rate_cats, bool invariant_sites,
          uint64_t seed, bool early_stop)
      : model_t(t, msa, rate_cats,
                std::vector<rate_category::rate_category_e>(
                    rate_cats.size(), rate_category::MEAN),
                invariant_sites, seed, early_stop){};
  ~model_t();
  double compute_lh(const root_location_t &root_location);
  double compute_lh_root(const root_location_t &root);
  dlh_t compute_dlh(const root_location_t &root_location);
  root_location_t optimize_alpha(const root_location_t &root, double atol);
  std::pair<root_location_t, double> optimize_root_location(size_t min_roots,
                                                            double root_ratio);

  std::pair<root_location_t, double> optimize_all(size_t min_roots,
                                                  double root_ratio,
                                                  double atol, double pgtol,
                                                  double brtol, double factor);

  std::pair<root_location_t, double>
  exhaustive_search(double atol, double pgtol, double brtol, double factor);

  const rooted_tree_t &rooted_tree(const root_location_t &root);
  const rooted_tree_t &virtual_rooted_tree(const root_location_t &root);
  const rooted_tree_t &unrooted_tree();

  void initialize_partitions(const std::vector<msa_t> &);
  void initialize_partitions_uniform_freqs(const std::vector<msa_t> &);
  std::string subst_string() const;

  std::vector<std::pair<root_location_t, double>> suggest_roots();
  std::vector<std::pair<root_location_t, double>> suggest_roots(size_t min,
                                                                double ratio);
  std::vector<root_location_t> suggest_roots_random(size_t min, double ratio);
  std::vector<double> compute_all_root_lh();
  void set_subst_rates(size_t, const model_params_t &);
  void set_freqs(size_t, const model_params_t &);
  void assign_indicies(const std::vector<size_t> &);
  void assign_indicies(size_t, size_t);
  void assign_indicies();
  void assign_indicies_by_rank(size_t, size_t);

private:
  std::pair<root_location_t, double>
  bisect(const root_location_t &beg, dlh_t d_beg, const root_location_t &end,
         dlh_t d_end, double atol, size_t depth);
  std::pair<root_location_t, double> brents(root_location_t beg, dlh_t d_beg,
                                            root_location_t end, dlh_t d_end,
                                            double atol);
  void set_subst_rates_random(size_t, const msa_t &);
  void set_subst_rates_random(size_t, size_t);
  void set_subst_rates_uniform();
  void set_gamma_weights(size_t, model_params_t);
  void set_gamma_rates(size_t);
  void set_gamma_rates(size_t, const model_params_t &);
  void set_gamma_rates_mean(size_t);
  void set_gamma_rates_mean(size_t, double);
  void set_gamma_rates_median(size_t);
  void set_gamma_rates_median(size_t, double);
  void set_gamma_rates_free(size_t);
  void set_gamma_rates_free(size_t, model_params_t);
  void update_invariant_sites(size_t);
  void set_tip_states(size_t, const msa_t &);
  void set_empirical_freqs(size_t);
  void set_empirical_freqs();
  void set_freqs_all_free(size_t, model_params_t);
  void move_root(const root_location_t &new_root);
  void update_pmatrices(const std::vector<unsigned int> &pmatrix_indices,
                        const std::vector<double> &branch_lengths);
  double bfgs_rates(model_params_t &initial_rates, const root_location_t &rl,
                    size_t partition_index, double pgtol, double factor);
  double bfgs_freqs(model_params_t &initial_rates, const root_location_t &rl,
                    size_t partition_index, double pgtol, double factor);
  double gd_rates(model_params_t &initial_rates, const root_location_t &rl,
                  size_t partition_index);
  double gd_freqs(model_params_t &initial_rates, const root_location_t &rl,
                  size_t partition_index);
  double bfgs_gamma_rates(model_params_t &intial_alpha,
                          const root_location_t &rl, size_t partition_index,
                          double pgtol, double factor);
  double gd_gamma_rates(model_params_t &intial_alpha, const root_location_t &rl,
                        size_t partition_index);
  double bfgs_gamma_weights(model_params_t &intial_alpha,
                            const root_location_t &rl, size_t partition_index,
                            double pgtol, double factor);
  double gd_gamma_weights(model_params_t &intial_alpha,
                          const root_location_t &rl, size_t partition_index);

  std::pair<size_t, size_t> compute_chunk_size_mod(size_t num_tasks) const;

#ifdef MPI_VERSION
  void gather_exhaustive_results(
      std::vector<std::pair<root_location_t, double>> &);
#endif

  rooted_tree_t _tree;
  std::vector<pll_partition_t *> _partitions;
  std::vector<rate_category::rate_category_e> _rate_category_types;
  std::vector<double> _partition_weights;
  std::vector<model_params_t> _rate_rates;
  std::vector<model_params_t> _rate_weights;
  std::vector<std::vector<unsigned int>> _param_indicies;
  std::vector<size_t> _assigned_idx;
  std::minstd_rand _random_engine;
  bool _invariant_sites;
  uint64_t _seed;
  bool _early_stop;
  /*
   * Only one submodel will be used for the time being. If there is desire for
   * more, we can add support for more models..
   */
  static constexpr unsigned int _submodels = 1;
};

#endif
