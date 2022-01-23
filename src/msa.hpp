#ifndef RD_MSA_HPP_
#define RD_MSA_HPP_

extern "C" {
#include <corax/corax.h>
}

#include "debug.h"
#include "util.hpp"
#include <string>
#include <unordered_set>
#include <vector>

typedef std::vector<partition_info_t> msa_partitions_t;

corax_msa_t     *parse_msa_file(const std::string &msa_filename);
msa_partitions_t parse_partition_file(const std::string &filename);
partition_info_t parse_partition_info(const std::string &line);
model_info_t     parse_model_info(const std::string &line);

class msa_t {
public:
  msa_t(const std::string   &msa_filename,
        const corax_state_t *map      = corax_map_nt,
        unsigned int         states   = 4,
        bool                 compress = true) :
      _msa{nullptr}, _map(map), _weights{nullptr}, _states(states) {
    _msa = parse_msa_file(msa_filename);
    if (compress) {
      _weights = corax_compress_site_patterns(
          _msa->sequence, map, _msa->count, &(_msa->length));
    } else {
      _weights = (unsigned int *)malloc(
          sizeof(unsigned int) * static_cast<unsigned long>(_msa->length));
      for (int i = 0; i < _msa->length; ++i) { _weights[i] = 1; }
    }
  };
  msa_t(const msa_t &other, const partition_info_t &part);
  msa_t(const msa_t &) = delete;
  msa_t(msa_t &&other) :
      _msa(other._msa),
      _map(other._map),
      _weights(other._weights),
      _states(other._states) {}

  char                *sequence(int) const;
  char                *label(int) const;
  unsigned int        *weights() const;
  unsigned int         total_weight() const;
  const corax_state_t *map() const;
  unsigned int         states() const;
  int                  count() const;
  unsigned int         length() const;
  std::vector<msa_t>   partition(const msa_partitions_t &) const;

  void compress();

  bool constiency_check(std::unordered_set<std::string>) const;
  void valid_data() const;

  ~msa_t();

private:
  corax_msa_t         *_msa;
  const corax_state_t *_map;
  unsigned int        *_weights;
  unsigned int         _states;
};

#endif
