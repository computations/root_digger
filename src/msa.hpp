#ifndef RD_MSA_HPP_
#define RD_MSA_HPP_

extern "C" {
#include <libpll/pll.h>
}

#include "debug.h"
#include <string>
#include <unordered_set>
#include <vector>

namespace param_type {
enum param_type_e { emperical, estimate, equal, user };
}

struct freq_opts_t {
  param_type::param_type_e type;
};

struct invar_opts_t {
  param_type::param_type_e type;
  float user_prop;
};

namespace rate_category {
enum rate_category_e { MEDIAN, MEAN, FREE };
}

struct ratehet_opts_t {
  param_type::param_type_e type;
  rate_category::rate_category_e rate_category_type;
  size_t rate_cats = 0;
  bool alpha_init = false;
  double alpha;
};

namespace asc_bias_type {
enum asc_bias_type_e { lewis, fels, stam };
}

struct asc_bias_opts_t {
  asc_bias_type::asc_bias_type_e type;
  double fels_weight;
  std::vector<double> stam_weights;
};

struct model_info_t {
  size_t states;
  std::string subst_str;
  freq_opts_t freq_opts;
  invar_opts_t invar_opts;
  ratehet_opts_t ratehet_opts;
  asc_bias_opts_t asc_opts;
};

struct partition_info_t {
  std::vector<std::pair<size_t, size_t>> parts;
  std::string model_name;
  std::string partition_name;
  model_info_t model;
};

/*
static const char *dna_models[] = {
    "JC",   "K80",    "F81",   "HKY",    "TN93ef", "TN93",   "K81",  "K81uf",
    "TPM2", "TPM2uf", "TPM3",  "TPM3uf", "TIM1",   "TIM1uf", "TIM2", "TIM2uf",
    "TIM3", "TIM3uf", "TVMef", "TVM",    "SYM",    "GTR"};

static const char *protein_models[] = {
    "Blosum62", "cpREV", "Dayhoff", "DCMut",       "DEN",    "FLU",
    "HIVb",     "HIVw",  "JTT",     "JTT - DCMut", "LG",     "mtART",
    "mtMAM",    "mtREV", "mtZOA",   "PMB",         "rtREV",  "stmtREV",
    "VT",       "WAG",   "LG4M ",   "LG4X ",       "PROTGTR"};

static const char *bin_models[] = {"BIN"};
*/

typedef std::vector<partition_info_t> msa_partitions_t;

pll_msa_t *parse_msa_file(const std::string &msa_filename);
msa_partitions_t parse_partition_file(const std::string &filename);
partition_info_t parse_partition_info(const std::string &line);
model_info_t parse_model_info(const std::string &line);

class msa_t {
public:
  msa_t(const std::string &msa_filename, const pll_state_t *map = pll_map_nt,
        unsigned int states = 4, bool compress = true)
      : _msa{nullptr}, _map(map), _weights{nullptr}, _states(states) {
    _msa = parse_msa_file(msa_filename);
    if (compress) {
      _weights = pll_compress_site_patterns(_msa->sequence, map, _msa->count,
                                            &(_msa->length));
    } else {
      _weights = (unsigned int *)malloc(
          sizeof(unsigned int) * static_cast<unsigned long>(_msa->length));
      for (int i = 0; i < _msa->length; ++i) {
        _weights[i] = 1;
      }
    }
  };
  msa_t(const msa_t &other, const partition_info_t &part);
  msa_t(const msa_t &) = delete;
  msa_t(msa_t &&other)
      : _msa(other._msa), _map(other._map), _weights(other._weights),
        _states(other._states) {}

  char *sequence(int) const;
  char *label(int) const;
  unsigned int *weights() const;
  unsigned int total_weight() const;
  const pll_state_t *map() const;
  unsigned int states() const;
  int count() const;
  unsigned int length() const;
  void compress();
  std::vector<msa_t> partition(const msa_partitions_t &) const;

  bool constiency_check(std::unordered_set<std::string>) const;
  void valid_data() const;

  ~msa_t();

private:
  pll_msa_t *_msa;
  const pll_state_t *_map;
  unsigned int *_weights;
  unsigned int _states;
};

#endif
