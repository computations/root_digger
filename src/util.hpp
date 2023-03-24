#ifndef RD_UTIL_H
#define RD_UTIL_H

#include <filesystem>
#include <random>

/* The following functions are taken from RAxML-NG. I did some of the initial
 * work, but the final implementation is from Alexey Kozlov
 */

#ifdef _OPENMP
bool sysutil_dir_exists(const std::string &dname);

std::string build_path(size_t cpu_number);

void get_cpuid(int32_t out[4], int32_t x);

size_t read_id_from_file(const std::string &filename);

size_t get_numa_node_id(const std::string &cpu_path);

size_t get_core_id(const std::string &cpu_path);

size_t get_physical_core_count(size_t n_cpu);

bool ht_enabled();

size_t sysutil_get_cpu_cores();
#else
size_t sysutil_get_cpu_cores();
#endif

std::string combine_argv_argc(int argv, char **argc);

typedef std::vector<double> model_params_t;

enum class param_type { emperical, estimate, equal, user };

struct freq_opts_t {
  param_type type;
};

struct invar_opts_t {
  param_type type;
  float      user_prop;
};

enum class rate_category { MEDIAN, MEAN, FREE };

struct ratehet_opts_t {
  ratehet_opts_t() = default;
  ratehet_opts_t(size_t rc) :
      type{param_type::estimate},
      rate_category_type{rate_category::MEAN},
      rate_cats{rc},
      alpha_init{false},
      alpha{1.0} {}

  param_type    type;
  rate_category rate_category_type;
  size_t        rate_cats  = 0;
  bool          alpha_init = false;
  double        alpha;

  bool operator==(const ratehet_opts_t &other) const {
    return type == other.type && rate_category_type == other.rate_category_type
           && rate_cats == other.rate_cats && alpha_init == other.alpha_init
           && alpha == other.alpha;
  }
};

enum class asc_bias_type { lewis, fels, stam };

enum class initial_root_strategy_t {
  random,
  midpoint,
  modified_mad,
};

struct asc_bias_opts_t {
  asc_bias_type       type;
  double              fels_weight;
  std::vector<double> stam_weights;
};

struct model_info_t {
  size_t          states;
  std::string     subst_str;
  freq_opts_t     freq_opts;
  invar_opts_t    invar_opts;
  ratehet_opts_t  ratehet_opts;
  asc_bias_opts_t asc_opts;
};

struct partition_info_t {
  std::vector<std::pair<size_t, size_t>> parts;
  std::string                            model_name;
  std::string                            partition_name;
  model_info_t                           model;
};

struct partition_parameters_t {
  model_params_t subst_rates;
  model_params_t freqs;
  model_params_t gamma_alpha;
  model_params_t gamma_weights;
};

struct rd_result_t {
  size_t root_id;
  double llh;
  double alpha;
};

class initialized_flag_t {
public:
  enum class initial_behavior {
    uninitalized,
    initialized_true,
    initialized_false,
  };
  initialized_flag_t() : value(initial_behavior::uninitalized){};

  initialized_flag_t(const initial_behavior &v) : value(v) {}
  /*
  initialized_flag_t &operator=(const initialized_flag_t &rhs) {
    value = rhs.value;
    return *this;
  }
  */

  bool operator==(const initialized_flag_t &rhs) const {
    return rhs.value == value;
  }
  bool operator!=(const initialized_flag_t &rhs) const {
    return rhs.value != value;
  }
  bool initalized() const { return value != initial_behavior::uninitalized; }

  bool convert_with_default(bool default_value) const {
    if (value == initial_behavior::uninitalized) return default_value;
    return value == initial_behavior::initialized_true;
  }

private:
  initial_behavior value;
};

struct cli_options_t {
  std::filesystem::path       msa_filename;
  std::filesystem::path       tree_filename;
  std::filesystem::path       prefix;
  std::filesystem::path       prefix_dir;
  std::string                 model_filename;
  std::string                 freqs_filename;
  std::filesystem::path       partition_filename;
  std::string                 data_type;
  std::string                 model_string;
  std::vector<ratehet_opts_t> rate_cats       = {1};
  uint64_t                    seed            = std::random_device()();
  size_t                      min_roots       = 1;
  size_t                      threads         = 0;
  double                      root_ratio      = 0.01;
  double                      abs_tolerance   = 1e-7;
  double                      factor          = 1e4;
  double                      br_tolerance    = 1e-12;
  double                      bfgs_tol        = 1e-7;
  unsigned int                states          = 4;
  bool                        silent          = false;
  bool                        exhaustive      = false;
  bool                        echo            = false;
  bool                        invariant_sites = false;
  bool                        clean           = false;
  initialized_flag_t          early_stop;

  initial_root_strategy_t initial_root_strategy = {
      initial_root_strategy_t::modified_mad};

  bool operator==(const cli_options_t &other) const {
    return msa_filename == other.msa_filename
           && tree_filename == other.tree_filename && prefix == other.prefix
           && prefix_dir == other.prefix_dir
           && model_filename == other.model_filename
           && freqs_filename == other.freqs_filename
           && partition_filename == other.partition_filename
           && data_type == other.data_type && model_string == other.model_string
           && rate_cats == other.rate_cats && seed == other.seed
           && threads == other.threads && root_ratio == other.root_ratio
           && abs_tolerance == other.abs_tolerance && factor == other.factor
           && br_tolerance == other.br_tolerance && bfgs_tol == other.bfgs_tol
           && states == other.states && exhaustive == other.exhaustive
           && echo == other.echo && invariant_sites == other.invariant_sites
           && early_stop == other.early_stop
           && initial_root_strategy == other.initial_root_strategy;
  }
  bool operator!=(const cli_options_t &other) const {
    return !(*this == other);
  }
};

#endif
