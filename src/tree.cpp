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
    _root_clv_index = 0;
    _root_scaler_index = 0;
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
    _root_clv_index = other._root_clv_index;
    _root_scaler_index = other._root_scaler_index;
  }
  return *this;
}

rooted_tree_t &rooted_tree_t::operator=(const rooted_tree_t &other) {
  if (_tree)
    pll_utree_destroy(_tree, nullptr);
  if (other._tree)
    _tree = pll_utree_clone(other._tree);
  else
    _tree = nullptr;
  generate_root_locations();
  return *this;
}

root_location_t rooted_tree_t::root_location(size_t index) const {
  if (index < _roots.size())
    return _roots[index];
  throw std::invalid_argument("Requested a root that does not exist");
}

size_t rooted_tree_t::root_count() const { return _roots.size(); }

std::vector<pll_unode_t *> rooted_tree_t::full_tree_traverse() const {
  return full_tree_traverse(_tree->vroot);
}

std::vector<pll_unode_t *>
rooted_tree_t::full_tree_traverse(pll_unode_t *vroot) const {
  std::vector<pll_unode_t *> trav_buf(_tree->edge_count);
  /*
   * pll_utree_traverse only puts each edge in once, so we can just use it to
   * produce a list of edges
   */
  auto edge_cb = [](pll_unode_t *node) -> int {
    if (node->next)
      return PLL_SUCCESS;
    return PLL_FAILURE;
  };

  unsigned int trav_size = 0;

  pll_utree_traverse(vroot, PLL_TREE_TRAVERSE_POSTORDER, edge_cb,
                     trav_buf.data(), &trav_size);
  trav_buf.resize(trav_size);
  return trav_buf;
}

void rooted_tree_t::generate_root_locations() {
  _roots.clear();
  if (_roots.capacity() < _tree->edge_count) {
    _roots.reserve(_tree->edge_count);
  }

  auto trav_buf = full_tree_traverse();

  for (size_t i = 0; i < trav_buf.size(); ++i) {
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

unsigned int rooted_tree_t::root_clv_index() const { return _root_clv_index; }
unsigned int rooted_tree_t::root_scaler_index() const {
  return _root_clv_index;
}

/*
 * Creates a map from tip label to clv index. Primarily used to match the
 * clv buffers with the corriesponding alignment.
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

std::tuple<std::vector<pll_operation_t>, std::vector<unsigned int>,
           std::vector<double>>
rooted_tree_t::generate_operations(const root_location_t &root) const {
  std::vector<pll_operation_t> ops;

  std::vector<pll_unode_t *> trav_buf = full_tree_traverse(root.edge);

  /* Reserve an aditional operation for the "extra" root clv */
  ops.resize(_tree->inner_count + 1);

  std::vector<double> branches(_tree->edge_count + 2);
  std::vector<unsigned int> pmatrix_indicies(_tree->edge_count + 2);

  unsigned int ops_size = 0;
  unsigned int matrix_count = 0;
  pll_utree_create_operations(trav_buf.data(), trav_buf.size(), branches.data(),
                              pmatrix_indicies.data(), ops.data(),
                              &matrix_count, &ops_size);

  pmatrix_indicies.resize(matrix_count + 2);

  /* manually construct the final operation */
  auto root_op_it = ops.end() - 1;

  root_op_it->parent_clv_index = _tree->inner_count + _tree->tip_count;
  root_op_it->parent_scaler_index = _tree->inner_count + _tree->tip_count;

  root_op_it->child1_clv_index = (root_op_it - 1)->parent_clv_index;
  root_op_it->child1_scaler_index = (root_op_it - 1)->parent_scaler_index;
  root_op_it->child1_matrix_index = _tree->edge_count;

  root_op_it->child2_clv_index = (root_op_it - 2)->parent_clv_index;
  root_op_it->child2_scaler_index = (root_op_it - 2)->parent_scaler_index;
  root_op_it->child2_matrix_index = _tree->edge_count + 1;

  auto root_branches_it = branches.end() - 1;
  *root_branches_it = root.brlen();
  *(root_branches_it - 1) = root.brlen_compliment();

  auto root_pmatrix_indicies_it = pmatrix_indicies.end() - 1;
  *root_pmatrix_indicies_it = _tree->edge_count + 1;
  *(root_pmatrix_indicies_it - 1) = _tree->edge_count;

  return std::make_tuple(ops, pmatrix_indicies, branches);
}
