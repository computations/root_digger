#ifndef __RD_MSA_HPP_
#define __RD_MSA_HPP_

extern "C" {
#include <libpll/pll.h>
}

#include <string>
#include <unordered_set>
#include <vector>

struct partition_info_t {
  size_t begin;
  size_t end;
  std::string model_name;
  std::string partition_name;
};

typedef std::vector<partition_info_t> msa_partitions_t;

pll_msa_t *parse_msa_file(const std::string &msa_filename);
msa_partitions_t parse_partition_file(const std::string &filename);
partition_info_t parse_partition_info(const std::string &line);

class msa_t {
public:
  msa_t(const std::string &msa_filename, const pll_state_t *map = pll_map_nt,
        unsigned int states = 4)
      : _msa{0}, _map(map), _weights{0}, _states(states) {
    _msa = parse_msa_file(msa_filename);
    _weights = pll_compress_site_patterns(_msa->sequence, map, _msa->count,
                                          &(_msa->length));
  };

  char *sequence(int) const;
  char *label(int) const;
  unsigned int *weights() const;
  const pll_state_t *map() const;
  unsigned int states() const;
  int count() const;
  int length() const;

  bool constiency_check(std::unordered_set<std::string>) const;

  ~msa_t();

private:
  pll_msa_t *_msa;
  const pll_state_t *_map;
  unsigned int *_weights;
  unsigned int _states;
};

#endif
