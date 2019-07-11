#include "debug.h"
#include "tree.hpp"
extern "C" {
#include <libpll/pll_tree.h>
}
#include <algorithm>
#include <iostream>
#include <unordered_set>

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
unsigned int rooted_tree_t::inner_count() const {
  return _tree->inner_count + 1;
}
unsigned int rooted_tree_t::branch_count() const {
  return _tree->tip_count * 2 - 2;
}

unsigned int rooted_tree_t::root_clv_index() const {
  return _tree->vroot->clv_index;
}
unsigned int rooted_tree_t::root_scaler_index() const {
  return _tree->vroot->scaler_index;
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
  debug_string("generating root locations");
  auto edges = full_traverse();

  std::unordered_set<pll_unode_t *> node_set;
  for (auto edge : edges) {
    if (node_set.find(edge) == node_set.end() &&
        node_set.find(edge->back) == node_set.end()) {
      node_set.insert(edge);
    }
  }

  _roots.reserve(node_set.size());
  for (auto edge : node_set) {
    _roots.push_back({edge, edge->length, 0.5});
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
  debug_print("rooting by node labeled: %s", root_location.label().c_str());
  if (root_location.edge == _tree->vroot) {
    update_root(root_location);
    return;
  }
  if (rooted()) {
    unroot();
  }
  /* make the roots */
  pll_unode_t *new_root_left = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
  pll_unode_t *new_root_right = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));

  new_root_left->next = new_root_right;
  new_root_right->next = new_root_left;

  pll_unode_t *left_child = root_location.edge;
  pll_unode_t *right_child = left_child->back;

  left_child->back = new_root_left;
  new_root_left->back = left_child;
  left_child->length = new_root_left->length = root_location.brlen();

  right_child->back = new_root_right;
  new_root_right->back = right_child;
  right_child->length = new_root_right->length =
      root_location.brlen_compliment();

  size_t new_size = _tree->inner_count + _tree->tip_count + 1;
  size_t total_unodes = _tree->inner_count * 3 + _tree->tip_count;

  _tree->nodes =
      (pll_unode_t **)realloc(_tree->nodes, sizeof(pll_unode_t *) * (new_size));

  _tree->nodes[new_size - 1] = new_root_left;
  _tree->inner_count += 1;
  _tree->edge_count += 1;
  _tree->vroot = new_root_left;

  new_root_left->clv_index = new_root_right->clv_index = new_size - 1;
  new_root_left->scaler_index = new_root_right->scaler_index =
      _tree->inner_count - 1;

  new_root_left->node_index = total_unodes + 1;
  new_root_left->pmatrix_index = left_child->pmatrix_index;

  new_root_right->node_index = total_unodes + 2;
  right_child->pmatrix_index = new_root_right->pmatrix_index =
      _tree->edge_count - 1;
  new_root_right->pmatrix_index = _tree->edge_count - 1;
}

void rooted_tree_t::update_root(root_location_t root) {
  if (root.edge != _tree->vroot) {
    throw std::runtime_error("Provided root doesn't match the current tree");
  }

  pll_unode_t *right_root = root.edge;
  pll_unode_t *left_root = root.edge->next;

  right_root->length = right_root->back->length = root.brlen();
  left_root->length = left_root->back->length = root.brlen_compliment();
}

void rooted_tree_t::unroot() {
  pll_unode_t *left_child = _tree->vroot->back;
  pll_unode_t *right_child = _tree->vroot->next->back;
  pll_unode_t *root_left = _tree->vroot;
  pll_unode_t *root_right = _tree->vroot->next;

  right_child->back = left_child;
  left_child->back = right_child;

  double new_length = root_left->length + root_right->length;
  right_child->length = left_child->length = new_length;

  root_left->back = nullptr;
  root_right->back = nullptr;

  root_left->next = nullptr;
  root_right->next = nullptr;

  free(root_left);
  free(root_right);

  _tree->vroot = left_child->next != nullptr ? left_child : right_child;
  if (_tree->vroot->next == nullptr) {
    throw std::runtime_error("unrooted to a tip");
  }

  size_t new_size = _tree->inner_count + _tree->tip_count - 1;
  _tree->nodes =
      (pll_unode_t **)realloc(_tree->nodes, sizeof(pll_unode_t *) * new_size);
  _tree->inner_count -= 1;
  _tree->edge_count -= 1;

  right_child->pmatrix_index = left_child->pmatrix_index;
}

bool rooted_tree_t::rooted() const {
  return _tree->vroot->next->next == _tree->vroot;
}

std::tuple<std::vector<pll_operation_t>, std::vector<unsigned int>,
           std::vector<double>>
rooted_tree_t::generate_operations(const root_location_t &new_root) {
  root_by(new_root);

  auto trav_buf = full_traverse();
  debug_string("traversal after root");
  for (auto node : trav_buf) {
    debug_print("traversal node label: %s, pmatrix index: %d, clv index: %d",
                (node->label != nullptr ? node->label : "null"),
                node->pmatrix_index, node->clv_index);
  }

  std::vector<pll_operation_t> ops(trav_buf.size());
  std::vector<unsigned int> pmatrix_indices(trav_buf.size());
  std::vector<double> branch_lengths(trav_buf.size());

  unsigned int op_count = 0;
  unsigned int matrix_count = 0;

  pll_utree_create_operations(trav_buf.data(), trav_buf.size() - 1,
                              branch_lengths.data(), pmatrix_indices.data(),
                              ops.data(), &matrix_count, &op_count);

  ops.resize(op_count + 1);
  pmatrix_indices.resize(matrix_count);
  branch_lengths.resize(matrix_count);

  auto root_op_it = ops.end() - 1;
  pll_unode_t *root_node = *(trav_buf.end() - 1);
  root_op_it->parent_clv_index = root_node->clv_index;
  root_op_it->parent_scaler_index = root_node->scaler_index;

  root_op_it->child1_clv_index = root_node->back->clv_index;
  root_op_it->child1_scaler_index = root_node->back->scaler_index;
  root_op_it->child1_matrix_index = root_node->back->pmatrix_index;

  root_op_it->child2_clv_index = root_node->next->back->clv_index;
  root_op_it->child2_scaler_index = root_node->next->back->scaler_index;
  root_op_it->child2_matrix_index = root_node->next->back->pmatrix_index;

  return std::make_tuple(ops, pmatrix_indices, branch_lengths);
}

std::tuple<pll_operation_t, std::vector<unsigned int>, std::vector<double>>
rooted_tree_t::generate_derivative_operations(const root_location_t &root) {
  root_by(root);

  pll_operation_t op;
  std::vector<unsigned int> pmatrix_indices(2);
  std::vector<double> branch_lengths(2);

  pll_unode_t *vroot = _tree->vroot;

  op.parent_clv_index = root_clv_index();
  op.parent_scaler_index = root_scaler_index();

  op.child1_clv_index = vroot->back->clv_index;
  op.child1_matrix_index = vroot->back->pmatrix_index;
  op.child1_scaler_index = vroot->back->scaler_index;
  pmatrix_indices[0] = vroot->back->pmatrix_index;
  branch_lengths[0] = vroot->back->length;

  op.child2_clv_index = vroot->next->back->clv_index;
  op.child2_matrix_index = vroot->next->back->pmatrix_index;
  op.child2_scaler_index = vroot->next->back->scaler_index;
  pmatrix_indices[1] = vroot->next->back->pmatrix_index;
  branch_lengths[1] = vroot->next->back->length;

  return std::make_tuple(op, pmatrix_indices, branch_lengths);
}

std::string rooted_tree_t::newick() const {
  char *newick_string = pll_utree_export_newick(_tree->vroot, nullptr);
  std::string ret{newick_string};
  free(newick_string);
  return ret;
}

void rooted_tree_t::show_tree() const {
  pll_utree_show_ascii(_tree->vroot,
                       PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_BRANCH_LENGTH);
}

bool rooted_tree_t::sanity_check() const {
  auto nodes = full_traverse();
  std::sort(nodes.begin(), nodes.end(), [](pll_unode_t *a, pll_unode_t *b) {
    return a->length < b->length;
  });
  double median = (nodes[nodes.size() / 2]->length +
                   nodes[(nodes.size() + 1) / 2]->length) /
                  2.0;
  if (median * 10 < (*(nodes.end() - 1))->length) {
    return false;
  }
  return true;
}
