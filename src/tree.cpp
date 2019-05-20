#include "tree.hpp"

pll_utree_t *parse_tree_file(const std::string &tree_filename) {
  return pll_utree_parse_newick(tree_filename.c_str());
}
