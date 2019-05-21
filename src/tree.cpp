#include "tree.hpp"

pll_utree_t *parse_tree_file(const std::string &tree_filename) {
  return pll_utree_parse_newick_unroot(tree_filename.c_str());
}

rooted_tree_t::rooted_tree_t(const rooted_tree_t &other) {
  if (other._tree) {
    _tree = pll_utree_clone(other._tree);
    if (other._roots.size() > 0) {
      generate_root_locations();
    }
  } else {
    _tree = nullptr;
  }
}

rooted_tree_t::~rooted_tree_t() {
  if (_tree)
    pll_utree_destroy(_tree, nullptr);
}

rooted_tree_t &rooted_tree_t::operator=(rooted_tree_t &&other) {
  if (this != &other) {
    if (_tree)
      pll_utree_destroy(_tree, nullptr);
    _tree = other._tree;
    other._tree = nullptr; // ensure that what is now _ours_ stays _ours_!
    _roots = std::move(other._roots);
  }
  return *this;
}

rooted_tree_t &rooted_tree_t::operator=(const rooted_tree_t &other) {
  if(_tree)
    pll_utree_destroy(_tree, nullptr);
  if(other._tree)
    _tree = pll_utree_clone(other._tree);
  else
    _tree = nullptr;
  _roots = other._roots;
  return *this;
}

void rooted_tree_t::generate_root_locations() {
  _roots.clear();
  if (_roots.capacity() < _tree->edge_count) {
    _roots.reserve(_tree->edge_count);
  }

  /*
   * pll_utree_traverse only puts each edge in once, so we can just use it to
   * produce a list of edges
   */
  auto edge_cb = [](pll_unode_t *node) -> int {
    if (node->next)
      return PLL_SUCCESS;
    return PLL_FAILURE;
  };

  pll_unode_t **trav_buf =
      (pll_unode_t **)malloc(sizeof(pll_unode_t *) * (_tree->edge_count));

  unsigned int trav_size = 0;

  pll_utree_traverse(_tree->vroot, PLL_TREE_TRAVERSE_POSTORDER, edge_cb,
                     trav_buf, &trav_size);

  for (size_t i = 0; i < trav_size; ++i) {
    _roots.push_back(root_location_t{trav_buf[i], 0.5});
  }
}

/*
 * We are going to cheat the partition, where we extra CLVS that don't have
 * corriesponding nodes on the tree. These will be used to fake a root, so we
 * don't need to keep inserting and removing nodes. So, the number of branches
 * (and clvs and pmatrices) is going to be 2 more than the typical unrooted tree
 * case.
 */
unsigned int rooted_tree_t::branches() const {
  return _tree->tip_count * 2 - 3 + 2;
}

unsigned int rooted_tree_t::inner_count() const { return _tree->inner_count; }

unsigned int rooted_tree_t::tip_count() const { return _tree->tip_count; }

/*
 * Creates a map from tip label to clv index. Primarily used to match the clv
 * buffers with the corriesponding alignment.
 */
std::unordered_map<std::string, unsigned int> rooted_tree_t::label_map() const {
  std::unordered_map<std::string, unsigned int> label_map;

  for (unsigned int i = 0; i < tip_count(); ++i) {
    pll_unode_t *cur_node = _tree->nodes[i];
    if (cur_node->label) {
      label_map[cur_node->label] = cur_node->clv_index;
    }
  }
  return label_map;
}
