#include "model.hpp"
#include <fstream>
#include <unordered_map>

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

// TODO: compress msa
model_t::model_t(const model_params_t &rate_parameters, pll_utree_t *tree,
                 pll_msa_t *msa, unsigned int states) {

  _tree = pll_utree_clone(tree);
  /*
   * We are going to cheat the partition, where we extra CLVS that don't have
   * corriesponding nodes on the tree. These will be used to fake a root, so we
   * don't need to keep inserting and removing nodes. So, the number of branches
   * (and clvs and pmatrices) is going to be 2 more than the typical unrooted
   * tree case.
   */
  unsigned int branches = _tree->tip_count * 2 - 3 + 2;

  /*
   * Only one submodel will be used for the time being.
   */
  constexpr unsigned int submodels = 1;

  /*
   * For now, nonrev only supports the "generic" cpu case. I might be able to
   * boost it up to the other simd architectures for CLVS, but that is untested
   * at this moment.
   */
  unsigned int attributes = 0;
  attributes |= PLL_ATTRIB_ARCH_CPU;
  // attributes |= PLL_ATTRIB_NONREV;

  _partition = pll_partition_create(_tree->tip_count, branches, states,
                                    msa->length, submodels, branches, submodels,
                                    branches, attributes);
  pll_set_subst_params(_partition, 0, rate_parameters.data());

  /* build a label map */
  auto full_trav_cb = [](pll_unode_t *) -> int { return PLL_SUCCESS; };

  pll_unode_t **trav_buffer = (pll_unode_t **)malloc(
      sizeof(pll_unode_t *) * (_tree->tip_count + _tree->inner_count));
  unsigned int trav_size = 0;
  pll_utree_traverse(_tree->vroot, PLL_TREE_TRAVERSE_POSTORDER, full_trav_cb,
                     trav_buffer, &trav_size);

  std::unordered_map<std::string, unsigned int> label_map;

  for (unsigned int i = 0; i < trav_size; ++i) {
    if (trav_buffer[i]->label) {
      label_map[trav_buffer[i]->label] = trav_buffer[i]->clv_index;
    }
  }

  /* use the label map to assign tip states in the partition */
  for (int i = 0; i < msa->count; ++i) {
    pll_set_tip_states(_partition, label_map.at(msa->label[i]), pll_map_nt,
                       msa->sequence[i]);
  }
}

model_t::~model_t() {
  pll_utree_destroy(_tree, nullptr);
  pll_partition_destroy(_partition);
}
