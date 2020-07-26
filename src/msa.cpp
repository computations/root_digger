#include "debug.h"
#include "msa.hpp"
#include "pll.h"
#include <cctype>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include <libpll/pll_msa.h>
}

typedef std::string::const_iterator string_iter_t;

pll_msa_t *parse_msa_file(const std::string &msa_filename) {

  debug_print(EMIT_LEVEL_DEBUG, "attempting to open msa: %s",
              msa_filename.c_str());
  if (pll_phylip_t *fd =
          pll_phylip_open(msa_filename.c_str(), pll_map_generic)) {
    pll_msa_t *pll_msa = nullptr;
    if ((pll_msa = pll_phylip_parse_interleaved(fd)) ||
        ((pll_phylip_rewind(fd) == PLL_SUCCESS) &&
         (pll_msa = pll_phylip_parse_sequential(fd)))) {
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

    if (sequences.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error(
          "The size of the sequence is too large to cast safely");
    }

    if (expected_sequence_length >
        static_cast<long int>(std::numeric_limits<int>::max())) {
      throw std::runtime_error(
          "The expected size of the sequence is too large to cast safely");
    }

    pll_msa->count = static_cast<int>(sequences.size());
    pll_msa->length = static_cast<int>(expected_sequence_length);
    if (pll_msa->count < 0) {
      throw std::runtime_error("The MSA had a negative count (overflow?)");
    }
    pll_msa->sequence = (char **)malloc(
        sizeof(char *) * static_cast<unsigned int>(pll_msa->count));
    pll_msa->label = (char **)malloc(sizeof(char *) *
                                     static_cast<unsigned int>(pll_msa->count));
    for (size_t i = 0; i < sequences.size(); ++i) {
      pll_msa->sequence[i] = sequences[i];
      pll_msa->label[i] = labels[i];
    }
    return pll_msa;
  }
  throw std::invalid_argument("Could not parse msa file");
}

/* Given a character, looks for the next instance of that character. If the
 * another character is encountered that isn't the specified one, the function
 * throws an exception. Ignores spaces The the returned iterator is one past the
 * character that was expected.
 */
static string_iter_t expect_next(string_iter_t itr, char c) {
  while (std::isspace(*itr)) {
    itr++;
  }
  if (std::tolower(c) != std::tolower(*itr)) {
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

static inline string_iter_t scan_word(string_iter_t iter) {
  auto start = iter;
  while (std::isalnum(*iter) || (*iter == '_') || (*iter == ':')) {
    ++iter;
  }
  if (iter == start) {
    throw std::runtime_error("Failed to find a word when scanning");
  }
  return iter;
}

static inline string_iter_t scan_float(string_iter_t iter) {
  auto start = iter;
  while (std::isdigit(*iter)) {
    debug_print(EMIT_LEVEL_DEBUG, "eating digit %c", *iter);
    iter++;
  }
  if (start == iter) {
    throw std::runtime_error("Expected float, but found something else");
  }

  /* now we expect either a '.', '+', or an 'e'. */
  while (true) {
    bool dot = false;
    switch (*iter) {
    case '.':
    case '+':
    case '-':
    case 'e':
      if (*iter == '.') {
        dot = true;
      }
      iter++;
      break;
    default:
      return iter;
    }
    auto cont_start = iter;
    while (*iter && std::isdigit(*++iter)) {
    }
    if (!dot && cont_start == iter) {
      throw std::runtime_error("Encountered a malformed floating point number");
    }
  }

  return iter;
}

static inline string_iter_t scan_integer(string_iter_t iter) {
  auto start = iter;
  while (*iter && std::isdigit(*++iter)) {
  }
  if (start == iter) {
    throw std::runtime_error("Expected a number, but found something else");
  }
  return iter;
}

static inline std::string parse_subst_str(string_iter_t &iter) {
  /* the matrix spec is anything before the first plus */
  auto end = scan_word(iter);
  std::string subst_string{iter, end};
  iter = end;
  return subst_string;
}

static inline string_iter_t
skip_to(string_iter_t iter, const std::function<bool(char)> &test_func) {
  while (!test_func(*++iter) && *iter) {
    debug_print(EMIT_LEVEL_DEBUG, "Skipping '%c'/%x", *iter, *iter);
    debug_print(EMIT_LEVEL_DEBUG, "test_func '%x", test_func(*iter));
  }
  debug_print(EMIT_LEVEL_DEBUG, "Finished skipping on '%c'/%x", *iter, *iter);
  if (!*iter) {
    throw std::runtime_error(
        std::string("Could not skip, found end of string intead"));
  }
  return iter;
}

static inline string_iter_t skip_space(string_iter_t iter) {
  while (std::isspace(*iter)) {
    ++iter;
  }
  return iter;
}

static inline freq_opts_t parse_freq_options(string_iter_t &iter) {
  freq_opts_t fi;
  iter++;
  switch (*iter) {
  default:
    fi.type = param_type::emperical;
    break;
  case 'c':
  case 'C':
    iter++;
    fi.type = param_type::emperical;
    break;
  case 'O':
  case 'o':
    fi.type = param_type::estimate;
    iter++;
    break;
  case 'E':
  case 'e':
    fi.type = param_type::equal;
    iter++;
    break;
  case 'U':
  case 'u':
    fi.type = param_type::user;
    /* Eat the input until we find the next '}' */
    iter = skip_to(iter, [](char d) -> bool { return d == '}'; });
    ++iter;
    break;
  }

  iter = skip_space(iter);
  return fi;
}

static inline invar_opts_t parse_invar_options(string_iter_t &iter) {
  invar_opts_t ii;
  iter++;
  switch (*iter) {
  default:
    ii.type = param_type::estimate;
    break;
  case 'O':
  case 'o':
    ii.type = param_type::estimate;
    iter++;
    break;
  case 'C':
  case 'c':
    ii.type = param_type::emperical;
    iter++;
    break;
  case 'U':
  case 'u':
    iter++;
    ii.type = param_type::user;
    iter = expect_next(iter, '{');
    auto end = scan_float(iter);
    debug_print(EMIT_LEVEL_DEBUG, "ending on %c/%x", *end, *end);
    ii.user_prop = std::stod(std::string{iter, end});
    iter = expect_next(end, '}');
    break;
  }

  iter = skip_space(iter);
  return ii;
}

static inline ratehet_opts_t parse_ratehet_options(string_iter_t &iter) {
  ratehet_opts_t ri;
  iter++;
  switch (*iter) {
  default:
    ri.type = param_type::estimate;
    ri.rate_category_type = rate_category::MEAN;
    if (std::isdigit(*iter)) {
      auto end = scan_integer(iter);
      int res = std::stoi(std::string(iter, end));
      debug_print(EMIT_LEVEL_DEBUG, "res.c_str(): %s",
                  std::string(iter, end).c_str());
      if (res < 0) {
        throw std::runtime_error(
            "Number of rate categories can not be less than zero");
      }
      ri.rate_cats = static_cast<size_t>(res);
      iter = end;
      debug_print(EMIT_LEVEL_DEBUG, "*iter: %c/%x", *iter, *iter);
      if (*iter == '{') {
        iter++;
        end = scan_float(iter);
        ri.alpha = stod(std::string(iter, end));
        ri.alpha_init = true;
        iter = expect_next(end, '}');
        ri.type = param_type::user;
        debug_print(EMIT_LEVEL_DEBUG, "*iter: %c/%x", *iter, *iter);
      }
    } else {
      ri.rate_cats = 4;
    }
    break;
  case 'A':
  case 'a':
    ri.rate_category_type = rate_category::MEDIAN;
    ri.type = param_type::estimate;
    ri.rate_cats = 4;
    iter++;
  }
  debug_print(EMIT_LEVEL_DEBUG, "*iter: %c/%x", *iter, *iter);
  iter = skip_space(iter);
  return ri;
}

static inline ratehet_opts_t parse_ratehet_free_options(string_iter_t &iter) {
  ratehet_opts_t rho;
  iter++;
  rho.type = param_type::estimate;
  rho.rate_category_type = rate_category::FREE;
  auto end = scan_integer(iter);

  int res = std::stoi(std::string(iter, end));
  if (res < 0) {
    throw std::runtime_error(
        "Number of rate categories can not be less than zero");
  }
  rho.rate_cats = static_cast<size_t>(res);
  iter = end;
  debug_print(EMIT_LEVEL_DEBUG, "iter: %c", *iter);
  if (*iter == '{') {
    debug_string(EMIT_LEVEL_WARNING, "Ignoring the user provided rate category "
                                     "weights as they are not supported");
    iter = skip_to(iter, [](char c) -> bool { return c == '}'; });
    iter = skip_to(iter, [](char c) -> bool { return c == '}'; });
    ++iter;
  }
  iter = skip_space(iter);
  return rho;
}

static inline asc_bias_opts_t parse_ascbias_options(string_iter_t &iter) {
  asc_bias_opts_t abo;
  iter = expect_next(iter, 'S');
  iter = expect_next(iter, 'C');
  iter = expect_next(iter, '_');
  switch (*iter) {
  case 'L':
  case 'l':
    abo.type = asc_bias_type::lewis;
    iter++;
    break;
  case 'F':
  case 'f': {
    abo.type = asc_bias_type::fels;
    iter = expect_next(iter, '{');
    auto end = scan_float(iter);
    abo.fels_weight = std::stod(std::string(iter, end));
    iter = expect_next(end, '}');
  } break;
  case 'S':
  case 's': {
    abo.type = asc_bias_type::stam;
    iter = expect_next(iter, '{');
    while (true) {
      auto end = scan_float(iter);
      abo.stam_weights.push_back(std::stod(std::string(iter, end)));
      if (*end == '}') {
        iter = end;
        iter++;
        break;
      }
      iter = expect_next(end, '/');
    }
  } break;
  }
  iter = skip_space(iter);
  return abo;
}

model_info_t parse_model_info(const std::string &model_string) {
  model_info_t model_info;
  /* parse the subst matrix
   * This is just a list of strings. We currently only support the GTR model and
   * binary data, so we need to mark that the model is not GTR or BIN, and then
   * warn that we are just going to ignore that information.
   */
  debug_print(EMIT_LEVEL_DEBUG, "model string: %s", model_string.c_str());
  auto iter = model_string.begin();
  auto end = scan_word(iter);
  model_info.subst_str = std::string(iter, end);
  debug_print(EMIT_LEVEL_DEBUG, "subst string: %s",
              model_info.subst_str.c_str());
  debug_print(EMIT_LEVEL_DEBUG, "end - iter: %lu", end - iter);
  iter = end;

  while (model_string.end() != iter && *iter) {
    debug_print(EMIT_LEVEL_DEBUG, "working on charcter: %c, %x", *iter, *iter);
    iter = expect_next(iter, '+');
    iter = skip_space(iter);
    switch (*iter) {
    case 'F':
    case 'f':
      model_info.freq_opts = parse_freq_options(iter);
      break;
    case 'I':
    case 'i':
      model_info.invar_opts = parse_invar_options(iter);
      break;
    case 'G':
    case 'g':
      model_info.ratehet_opts = parse_ratehet_options(iter);
      break;
    case 'A':
    case 'a':
      model_info.asc_opts = parse_ascbias_options(iter);
      break;
    case 'R':
    case 'r':
      model_info.ratehet_opts = parse_ratehet_free_options(iter);
      break;
    case 'M':
    case 'm':
      debug_string(EMIT_LEVEL_WARNING,
                   "The 'M' option in the model string is not supported by "
                   "root digger, and will be ignored.");
      iter++;
      break;
    }
  }
  return model_info;
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
  while (std::isalnum(*itr) || *itr == '+' || *itr == '{' || *itr == '}' ||
         *itr == '/' || *itr == '.') {
    itr++;
  }
  pi.model_name = std::string(start, itr);

  pi.model = parse_model_info(pi.model_name);
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

  do {
    if (*itr == ',')
      ++itr;
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
      auto tmp = std::stol(std::string(start, itr));
      if (tmp < 0) {
        throw std::runtime_error(
            "There was a negative number in the partition specification");
      }
      begin = static_cast<size_t>(tmp);
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
      auto tmp = std::stol(std::string(start, itr));
      if (tmp < 0) {
        throw std::runtime_error(
            "There was a negative number in the partition specification");
      }
      end = static_cast<size_t>(tmp);
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
  } while (*itr == ',');

  if (pi.model_name.empty()) {
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

msa_t::msa_t(const msa_t &other, const partition_info_t &partition) {
  _msa = (pll_msa_t *)malloc(sizeof(pll_msa_t));
  _msa->count = other.count();
  _msa->length = 0;
  for (auto range : partition.parts) {
    /*
     * Since the range specification is [first, second], we have to add one to
     * include the endpoint
     */
    size_t cur_partition_length = (range.second - range.first) + 1;
    if (cur_partition_length >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::runtime_error("Partition range is too large to cast safely");
    }
    _msa->length += static_cast<int>(cur_partition_length);
    debug_print(EMIT_LEVEL_DEBUG, "%d", _msa->length);
  }

  if (other.count() < 0) {
    throw std::runtime_error("MSA from which we are copy constructing has a "
                             "negative count (overflow?)");
  }

  _msa->sequence = (char **)malloc(static_cast<size_t>(sizeof(char *)) *
                                   static_cast<size_t>(other.count()));

  for (int i = 0; i < other.count(); ++i) {
    if (_msa->length <= 0) {
      throw std::runtime_error(
          "The size of the MSA is less than 0 (overflow?)");
    }

    // We need to make space for the terminating zero
    _msa->sequence[i] = (char *)malloc(static_cast<size_t>(sizeof(char)) *
                                       static_cast<size_t>(_msa->length + 1));
    size_t cur_idx = 0;
    char *other_sequence = other.sequence(i);
    for (auto range : partition.parts) {
      if (range.first == 0) {
        throw std::runtime_error(
            "Partition ranges start at 1, but we encountered a 0");
      }
      for (size_t j = range.first - 1; j < range.second; ++j) {
        _msa->sequence[i][cur_idx++] = other_sequence[j];
      }
    }
  }

  if (_msa->count < 0) {
    throw std::runtime_error(
        "The current MSA count is less than 0 (overflow?)");
  }
  _msa->label =
      (char **)calloc(sizeof(char *), static_cast<size_t>(_msa->count));
  for (int i = 0; i < _msa->count; ++i) {
    size_t label_size = strlen(other.label(i)) + 1;
    _msa->label[i] = (char *)malloc(sizeof(char) * label_size);
    strncpy(_msa->label[i], other.label(i), label_size);
  }
  _map = other.map();
  _states = other.states();
  _weights = nullptr;
  compress();
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

unsigned int msa_t::total_weight() const {
  if (!_weights) {
    return static_cast<unsigned int>(_msa->length);
  }
  unsigned int total = 0;
  for (int i = 0; i < _msa->length; ++i) {
    total += _weights[i];
  }
  return total;
}

const pll_state_t *msa_t::map() const { return _map; }

unsigned int msa_t::states() const { return _states; }

int msa_t::count() const { return _msa->count; }

unsigned int msa_t::length() const {
  return static_cast<unsigned int>(_msa->length);
}

void msa_t::compress() {
  if (_weights != nullptr) {
    free(_weights);
  }
  int new_length = _msa->length;
  _weights = pll_compress_site_patterns(_msa->sequence, _map, _msa->count,
                                        &new_length);
  if (!_weights) {
    throw std::runtime_error(std::string("PLL ERR: ") +
                             std::to_string(pll_errno) + " " + pll_errmsg);
  }
  _msa->length = new_length;
}

std::vector<msa_t> msa_t::partition(const msa_partitions_t &ps) const {
  std::vector<msa_t> parted_msa;
  parted_msa.reserve(ps.size());
  for (auto p : ps) {
    parted_msa.emplace_back(*this, p);
  }
  return parted_msa;
}

bool msa_t::constiency_check(std::unordered_set<std::string> labels) const {
  std::unordered_set<std::string> taxa;
  for (int i = 0; i < _msa->count; ++i) {
    taxa.insert(_msa->label[i]);
  }

  bool ret = true;

  // labels subset taxa
  for (const std::string &k : labels) {
    if (taxa.find(k) == taxa.end()) {
      debug_print(EMIT_LEVEL_ERROR, "Taxa %s in msa is not present on the tree",
                  k.c_str());
      ret = false;
    }
  }

  // taxa subset labels
  for (const std::string &k : taxa) {
    if (labels.find(k) == labels.end()) {
      debug_print(EMIT_LEVEL_ERROR, "Taxa %s on tree is not present in the msa",
                  k.c_str());
      ret = false;
    }
  }
  return ret;
}

void msa_t::valid_data() const {
  for (int i = 0; i < _msa->count; ++i) {
    for (int j = 0; j < _msa->length; ++j) {
      char c = _msa->sequence[i][j];
      if (c < 0) {
        throw std::runtime_error(
            std::string("Encountered an invalid character in sequence ") +
            std::to_string(i) + " at position " + std::to_string(j) + ".");
      }
      size_t idx = static_cast<size_t>(c);
      if (_map[idx] == 0) {
        throw std::runtime_error(
            std::string("Found unrecognized character sequence ") +
            std::to_string(i) + " position " + std::to_string(j) + ".");
      }
    }
  }
}

msa_t::~msa_t() {
  if (_msa)
    pll_msa_destroy(_msa);
  if (_weights)
    free(_weights);
}
