#include "msa.hpp"
#include <vector>
#include <stdexcept>

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
