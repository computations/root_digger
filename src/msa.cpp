#include "debug.h"
#include "msa.hpp"
#include <cctype>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

extern "C" {
#include <libpll/pll_msa.h>
}

pll_msa_t *parse_msa_file(const std::string &msa_filename) {

  if (pll_phylip_t *fd =
          pll_phylip_open(msa_filename.c_str(), pll_map_generic)) {
    pll_msa_t *pll_msa = nullptr;
    if ((pll_msa = pll_phylip_parse_interleaved(fd)) ||
        (pll_msa = pll_phylip_parse_sequential(fd))) {
      pll_phylip_close(fd);
      return pll_msa;
    } else {
      pll_phylip_close(fd);
    }
  }

  if (pll_fasta_t *fd = pll_fasta_open(msa_filename.c_str(), pll_map_fasta)) {
    char *label = nullptr;
    char *sequence = nullptr;
    long sequence_len = 0;
    long header_len = 0;
    long sequence_number = 0;
    long expected_sequence_length = 0;

    std::vector<char *> sequences;
    std::vector<char *> labels;

    while (pll_fasta_getnext(fd, &label, &header_len, &sequence, &sequence_len,
                             &sequence_number)) {
      if (expected_sequence_length == 0) {
        expected_sequence_length = sequence_len;
      } else if (expected_sequence_length != sequence_len) {
        throw std::invalid_argument("Sequences don't match in size");
      }
      sequences.push_back(sequence);
      labels.push_back(label);
    }
    pll_fasta_close(fd);
    pll_msa_t *pll_msa = (pll_msa_t *)malloc(sizeof(pll_msa_t));
    pll_msa->count = sequences.size();
    pll_msa->length = expected_sequence_length;
    pll_msa->sequence = (char **)malloc(sizeof(char *) * pll_msa->count);
    pll_msa->label = (char **)malloc(sizeof(char *) * pll_msa->count);
    for (size_t i = 0; i < sequences.size(); ++i) {
      pll_msa->sequence[i] = sequences[i];
      pll_msa->label[i] = labels[i];
    }
    return pll_msa;
  }
  throw std::invalid_argument("Could not parse msa file");
}

std::string::const_iterator expect_next(std::string::const_iterator itr,
                                        char c) {
  while (std::isspace(*itr)) {
    itr++;
  }
  if (c != *itr) {
    throw std::runtime_error(
        std::string("Failed to parse partition file, expected '") + c +
        "' got '" + *itr + "' instead");
  }
  ++itr;
  while (std::isspace(*itr)) {
    itr++;
  }
  return itr;
}

partition_info_t parse_partition_info(const std::string &line) {
  partition_info_t pi;
  auto itr = line.begin();

  /* skip spaces */
  while (std::isspace(*itr)) {
    itr++;
  }

  /* parse model name */
  auto start = itr;
  while (std::isalnum(*itr)) {
    itr++;
  }
  pi.model_name = std::string(start, itr);

  /* check for comma */
  itr = expect_next(itr, ',');

  /* parse partition name */
  start = itr;
  while (std::isalnum(*itr) || *itr == '_') {
    itr++;
  }
  pi.partition_name = std::string(start, itr);

  /* check for = */
  itr = expect_next(itr, '=');

  do{
    if(*itr == ',') ++itr;
    while (std::isspace(*itr)) {
      itr++;
    }

    /* parse begin */
    size_t begin = 0;
    size_t end = 0;
    start = itr;
    while (std::isdigit(*itr)) {
      itr++;
    }
    try {
      begin = std::stol(std::string(start, itr));
    } catch (...) {
      throw std::runtime_error(
          std::string("Failed to parse beginning partition number") +
          std::string(start, itr) + " start: '" + *start + "', itr: '" + *itr +
          "'");
    }

    /* check for - */
    itr = expect_next(itr, '-');

    /* parse begin */
    start = itr;
    while (std::isdigit(*itr)) {
      itr++;
    }
    try {
      end = std::stol(std::string(start, itr));
    } catch (...) {
      throw std::runtime_error(
          std::string("Failed to parse beginning partition number") +
          std::string(start, itr) + " start: '" + *start + "', itr: '" + *itr +
          "'");
    }
    if (end < begin) {
      throw std::runtime_error(std::string("The end index of the partition '") +
                               pi.partition_name +
                               "' comes before the beginning");
    }
    pi.parts.emplace_back(begin, end);
    while (std::isspace(*itr)) {
      itr++;
    }
  } while(*itr == ',');

  if(pi.model_name.empty()){
    throw std::runtime_error("Error, partition is missing a model name");
  }

  return pi;
}

/* Grammer:
 * <MODEL_NAME> , <PARTITION_NAME> = <BEGIN> - <END>[, <BEGIN> - <END>]*
 */
msa_partitions_t parse_partition_file(const std::string &filename) {
  std::ifstream partition_file{filename};
  msa_partitions_t parts;
  for (std::string line; std::getline(partition_file, line);) {
    parts.push_back(parse_partition_info(line));
  }
  return parts;
}

void msa_t::partition(const msa_partitions_t& parts){
}

char *msa_t::sequence(int index) const {
  if (index < _msa->count && !(index < 0))
    return _msa->sequence[index];
  throw std::out_of_range("Requested sequence does not exist");
}

char *msa_t::label(int index) const {
  if (index < _msa->count && !(index < 0))
    return _msa->label[index];
  throw std::out_of_range("Requested label does not exist");
}

unsigned int *msa_t::weights() const {
  if (_weights)
    return _weights;
  throw std::runtime_error("msa_t has no weights");
}

const pll_state_t *msa_t::map() const { return _map; }

unsigned int msa_t::states() const { return _states; }

int msa_t::count() const { return _msa->count; }

int msa_t::length() const { return _msa->length; }

bool msa_t::constiency_check(std::unordered_set<std::string> labels) const {
  std::unordered_set<std::string> taxa;
  for (int i = 0; i < _msa->count; ++i) {
    taxa.insert(_msa->label[i]);
  }

  bool ret = true;

  // labels subset taxa
  for (const std::string &k : labels) {
    if (taxa.find(k) == taxa.end()) {
      debug_print("Taxa %s in msa is not present on the tree", k.c_str());
      ret = false;
    }
  }

  // taxa subset labels
  for (const std::string &k : taxa) {
    if (labels.find(k) == labels.end()) {
      debug_print("Taxa %s on tree is not present in the msa", k.c_str());
      ret = false;
    }
  }
  return ret;
}

msa_t::~msa_t() {
  if (_msa)
    pll_msa_destroy(_msa);
  if (_weights)
    free(_weights);
}
