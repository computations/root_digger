#ifndef __RD_MSA_HPP_
#define __RD_MSA_HPP_

extern "C" {
#include <libpll/pll.h>
}

#include <string>
#include <vector>
#include <unordered_set>

pll_msa_t *parse_msa_file(const std::string &msa_filename);

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
