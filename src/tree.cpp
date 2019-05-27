#include "tree.hpp"

pll_utree_t *parse_tree_file(const std::string &tree_filename) {
  return pll_utree_parse_newick_unroot(tree_filename.c_str());
}

rooted_tree_t::rooted_tree_t(const rooted_tree_t &other) {
  _tree = pll_utree_clone(other._tree);
  generate_root_locations();
}

rooted_tree_t::~rooted_tree_t() {
  if (_tree != nullptr) {
    pll_utree_destroy(_tree, nullptr);
  }
}

rooted_tree_t &rooted_tree_t::operator=(rooted_tree_t &&other) {
  _tree = std::move(other._tree);
  _roots = std::move(other._roots);
  other._tree = nullptr;
  return *this;
}

rooted_tree_t &rooted_tree_t::operator=(const rooted_tree_t &other) {
  _tree = pll_utree_clone(other._tree);
  generate_root_locations();
  return *this;
}

root_location_t rooted_tree_t::root_location(size_t index) const {
  if (index > _roots.size()) {
    throw std::invalid_argument("Invalid index for roots on this tree");
  }
  return _roots[index];
}

size_t rooted_tree_t::root_count() const { return _roots.size(); }

unsigned int rooted_tree_t::tip_count() const { return _tree->tip_count; }
unsigned int rooted_tree_t::inner_count() const { return _tree->inner_count; }
unsigned int rooted_tree_t::branch_count() const {
  return _tree->tip_count * 2 - 2;
}

std::unordered_map<std::string, unsigned int> rooted_tree_t::label_map() const {
  std::unordered_map<std::string, unsigned int> label_map;
  for (unsigned int i = 0; i < tip_count(); ++i) {
    auto cur_node = _tree->nodes[i];
    label_map[cur_node->label] = cur_node->clv_index;
  }
  return label_map;
}

void rooted_tree_t::generate_root_locations() {
  auto edges = edge_traverse();
  _roots.resize(edges.size());

  for (size_t i = 0; i < edges.size(); ++i) {
    _roots[i] = {edges[i], 0.5};
  }
}

std::vector<pll_unode_t *> rooted_tree_t::edge_traverse() const {
  std::vector<pll_unode_t *> trav_buf(inner_count());
  unsigned int trav_size = 0;

  auto edge_trav_cb = [](pll_unode_t *node) -> int {
    if (node->next != nullptr) {
      return PLL_SUCCESS;
    }
    return PLL_FAILURE;
  };

  pll_utree_traverse(_tree->vroot, PLL_TREE_TRAVERSE_POSTORDER, edge_trav_cb,
                     trav_buf.data(), &trav_size);
  trav_buf.resize(trav_size);
  return trav_buf;
}

std::vector<pll_unode_t *> rooted_tree_t::full_traverse() const {
  std::vector<pll_unode_t *> trav_buf(tip_count() + inner_count());
  unsigned int trav_size = 0;

  auto full_trav_cb = [](pll_unode_t *) -> int { return PLL_SUCCESS; };

  pll_utree_traverse(_tree->vroot, PLL_TREE_TRAVERSE_POSTORDER, full_trav_cb,
                     trav_buf.data(), &trav_size);
  trav_buf.resize(trav_size);
  return trav_buf;
}

void rooted_tree_t::root_by(const root_location_t &root_location) {
  pll_unode_t *root_node_left = (pll_unode_t *)malloc(sizeof(pll_unode_t));
  pll_unode_t *root_node_right = (pll_unode_t *)malloc(sizeof(pll_unode_t));

  /* bind them as a "node" */
  root_node_left->next = root_node_right;
  root_node_right->next = root_node_left;

  /*insert them into the tree */
  pll_unode_t *old_back = root_location.edge->back;

  root_location.edge->back = root_node_left;
  root_node_left->back = root_location.edge;

  old_back->back = root_node_right;
  root_node_right->back = old_back;

  /* update the lengths */
  double left_length = root_location.brlen();
  double right_length = root_location.brlen_compliment();

  root_location.edge->length = root_node_left->length = left_length;
  old_back->length = root_node_right->length = right_length;

  /* update the tree */
  unsigned int saved_tip_count = tip_count();
  unsigned int saved_saved_inner = inner_count() + 1;
  free(_tree->nodes);
  free(_tree);
  _tree = pll_utree_wraptree_multi(root_node_left, saved_tip_count,
                                   saved_saved_inner);
  if (_tree == PLL_FAILURE) {
    throw std::runtime_error("Failed to wrap the tree after rooting");
  }
  pll_utree_reset_template_indices(_tree->vroot, saved_tip_count);
}

void rooted_tree_t::unroot() {
  /* check if the tree is already a binary tree */
  if (_tree->vroot != _tree->vroot->next->next) {
    return;
  }
  pll_unode_t *old_root_left = _tree->vroot;
  pll_unode_t *old_root_right = _tree->vroot -> next;

  pll_unode_t *left_child = old_root_left->back;
  pll_unode_t *right_child = old_root_right->back;

  double old_length = old_root_left->length + old_root_right->length;

  left_child->length = right_child->length = old_length;

  left_child->back = right_child;
  right_child->back = left_child;

  unsigned int saved_tip_count = tip_count();
  unsigned int saved_saved_inner = inner_count() - 1;
  free(_tree->nodes);
  free(_tree);
  _tree = pll_utree_wraptree_multi(left_child, saved_tip_count,
                                   saved_saved_inner);
  if (_tree == PLL_FAILURE) {
    throw std::runtime_error("Failed to wrap the tree after rooting");
  }
  pll_utree_reset_template_indices(_tree->vroot, saved_tip_count);

  free(old_root_right);
  free(old_root_left);
}

std::tuple<std::vector<pll_operation_t>, std::vector<unsigned int>,
           std::vector<double>>
rooted_tree_t::generate_operations(const root_location_t &new_root) {

  root_by(new_root);

  auto trav_buf = full_traverse();

  std::vector<pll_operation_t> ops(trav_buf.size());
  std::vector<unsigned int> pmatrix_indices(trav_buf.size());
  std::vector<double> branch_lengths(trav_buf.size());

  unsigned int op_count = 0;
  unsigned int matrix_count = 0;

  pll_utree_create_operations(trav_buf.data(), trav_buf.size(),
                              branch_lengths.data(), pmatrix_indices.data(),
                              ops.data(), &matrix_count, &op_count);

  unroot();

  return std::make_tuple(ops, pmatrix_indices, branch_lengths);
}
