#include "tree.hpp"
#include <stdexcept>
extern "C" {
#include <libpll/pll_tree.h>
}
#include <algorithm>
#include <limits>
#include <numeric>
#include <unordered_set>

pll_utree_t *parse_tree_file(const std::string &tree_filename) {
  return pll_utree_parse_newick_unroot(tree_filename.c_str());
}

static void clean_root_unode(pll_unode_t *node) {
  assert(!node->data);

  node->length = -1;
  // We _want_ this to fail if we ever use it after without setting it. So we
  // set it to a ridiculous value
  node->node_index = std::numeric_limits<unsigned int>::max();

  node->back = nullptr;
}

rooted_tree_t::rooted_tree_t(const rooted_tree_t &other) {
  if (other.rooted()) {
    throw std::runtime_error{"Attempted to copy a tree that is rooted"};
  }
  _tree   = pll_utree_clone(other._tree);
  _rooted = other._rooted;

  copy_root_locations(other);
  copy_annotations(other);

  if (!_rooted) { add_root_space(); }
}

rooted_tree_t::~rooted_tree_t() {
  auto deallocate_data = [](void *n) {
    assert(n != (void *)0xdeadbeef);
    if (n != nullptr) free(n);
  };
  if (_tree != nullptr) {
    // pll_utree_destroy works by traversing the tree. We just need to make sure
    // that the extra roots are actually part of the tree in order to not leak
    // memory.
    root_by(_roots[0]);
    pll_utree_destroy(_tree, deallocate_data);
  }
}

rooted_tree_t &rooted_tree_t::operator=(rooted_tree_t &&other) {
  _tree             = std::move(other._tree);
  _roots            = std::move(other._roots);
  _root_annotations = std::move(other._root_annotations);
  _rooted           = other._rooted;
  other._tree       = nullptr;
  return *this;
}

rooted_tree_t &rooted_tree_t::operator=(const rooted_tree_t &other) {
  if (other.rooted()) {
    throw std::runtime_error{"Attempted to copy a tree that is rooted"};
  }
  _tree   = pll_utree_clone(other._tree);
  _rooted = other._rooted;

  copy_root_locations(other);
  copy_annotations(other);

  if (!_rooted) { add_root_space(); }
  return *this;
}

root_location_t rooted_tree_t::root_location(size_t index) const {
  if (index > _roots.size()) {
    throw std::invalid_argument(
        std::string("Invalid index for roots on this tree: ")
        + std::to_string(index));
  }
  return _roots[index];
}

root_location_t rooted_tree_t::root_location(const std::string &name) const {
  for (auto rl : _roots) {
    if (rl.edge->label && name == rl.edge->label) { return rl; }
  }
  throw std::runtime_error{
      std::string{"Can't find the root location with label: "} + name};
}

root_location_t rooted_tree_t::root_location() const { return _current_rl; }

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
int rooted_tree_t::root_scaler_index() const {
  return _tree->vroot->scaler_index;
}

std::unordered_map<std::string, unsigned int> rooted_tree_t::label_map() const {
  std::unordered_map<std::string, unsigned int> label_map;
  if (_tree == nullptr) { return label_map; }
  for (unsigned int i = 0; i < tip_count(); ++i) {
    auto cur_node              = _tree->nodes[i];
    label_map[cur_node->label] = cur_node->clv_index;
  }
  return label_map;
}

std::unordered_set<std::string> rooted_tree_t::label_set() const {
  std::unordered_set<std::string> label_set;
  if (_tree == nullptr) { return label_set; }
  for (unsigned int i = 0; i < tip_count(); ++i) {
    label_set.insert(_tree->nodes[i]->label);
  }
  return label_set;
}

void rooted_tree_t::copy_root_locations(const rooted_tree_t &other) {
  _roots.clear();
  std::unordered_map<pll_unode_t *, size_t> id_map;
  id_map.reserve(other.root_count() * 2);
  for (const auto &r : other.roots()) {
    debug_print(EMIT_LEVEL_DEBUG, "mapping %p to %lu", (void *)r.edge, r.id);
    id_map[r.edge]       = r.id;
    id_map[r.edge->back] = r.id;
  }

  auto other_trav = other.full_traverse();
  auto our_trav   = full_traverse();

  if (other_trav.size() != our_trav.size()) {
    throw std::runtime_error("Traversal sizes didn't match during copy "
                             "constructor, something is seriously wrong");
  }

  std::unordered_set<size_t> consumed_ids;
  for (size_t i = 0; i < other_trav.size(); ++i) {
    auto id_map_res = id_map.find(other_trav[i]);
    if (id_map_res != id_map.end()
        && consumed_ids.find(id_map_res->second) == consumed_ids.end()) {
      auto edge = our_trav[i];
      _roots.push_back({edge, id_map_res->second, edge->length, 0.5});
      consumed_ids.insert(id_map_res->second);
    }
  }

  if (_roots.size() != other._roots.size()) {
    throw std::runtime_error{"We got the wrong number of roots after copy"};
  }
  sort_root_locations();
}

void rooted_tree_t::sort_root_locations() {
  std::sort(_roots.begin(),
            _roots.end(),
            [](const root_location_t &a, const root_location_t &b) {
              return a.id < b.id;
            });
}

void rooted_tree_t::generate_root_locations() {
  debug_string(EMIT_LEVEL_DEBUG, "generating root locations");
  auto edges = full_traverse();

  std::unordered_set<pll_unode_t *> node_set;
  _roots.reserve(_tree->inner_count + _tree->tip_count);
  node_set.reserve(_tree->inner_count + _tree->tip_count);
  size_t id = 0;
  for (auto &edge : edges) {
    if (node_set.find(edge) == node_set.end()
        && node_set.find(edge->back) == node_set.end()) {
      node_set.insert(edge);
      _roots.push_back({edge, id++, edge->length, 0.5});
    }
  }
}

void rooted_tree_t::add_root_space() {
  unsigned int new_size     = _tree->inner_count + _tree->tip_count + 1;
  unsigned int total_unodes = _tree->inner_count * 3 + _tree->tip_count;
  _tree->nodes =
      (pll_unode_t **)realloc(_tree->nodes, sizeof(pll_unode_t *) * new_size);
  assert(new_size > 0);
  pll_unode_t *new_root_left  = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
  pll_unode_t *new_root_right = (pll_unode_t *)calloc(1, sizeof(pll_unode_t));
  new_root_left->next         = new_root_right;
  new_root_right->next        = new_root_left;

  new_root_left->clv_index = new_root_right->clv_index = new_size;
  new_root_left->scaler_index = new_root_right->scaler_index =
      static_cast<int>(_tree->inner_count - 1);

  new_root_left->node_index = total_unodes + 1;

  new_root_right->node_index    = total_unodes + 2;
  new_root_right->pmatrix_index = _tree->edge_count - 1;

  _tree->nodes[new_size - 1] = new_root_left;
}

std::vector<pll_unode_t *> rooted_tree_t::edge_traverse() const {
  std::vector<pll_unode_t *> trav_buf(inner_count());
  unsigned int               trav_size = 0;

  auto edge_trav_cb = [](pll_unode_t *node) -> int {
    if (node->next != nullptr) { return PLL_SUCCESS; }
    return PLL_FAILURE;
  };

  pll_utree_traverse(_tree->vroot,
                     PLL_TREE_TRAVERSE_POSTORDER,
                     edge_trav_cb,
                     trav_buf.data(),
                     &trav_size);
  trav_buf.resize(trav_size);
  return trav_buf;
}

std::vector<pll_unode_t *> rooted_tree_t::full_traverse() const {
  std::vector<pll_unode_t *> trav_buf(tip_count() + inner_count());
  unsigned int               trav_size = 0;

  auto full_trav_cb = [](pll_unode_t *) -> int { return PLL_SUCCESS; };

  pll_utree_traverse(_tree->vroot,
                     PLL_TREE_TRAVERSE_POSTORDER,
                     full_trav_cb,
                     trav_buf.data(),
                     &trav_size);
  trav_buf.resize(trav_size);
  return trav_buf;
}

void rooted_tree_t::root_by(unsigned int root_id) { root_by(_roots[root_id]); }

void rooted_tree_t::root_by(const root_location_t &root_location) {
  debug_print(EMIT_LEVEL_DEBUG,
              "rooting by node labeled: %s",
              root_location.label().c_str());
  if (root_location.edge == _tree->vroot) {
    update_root(root_location);
    return;
  }
  if (rooted()) { unroot(); }
  unsigned int tree_size      = _tree->inner_count + _tree->tip_count + 1;
  pll_unode_t *new_root_left  = _tree->nodes[tree_size - 1];
  pll_unode_t *new_root_right = new_root_left->next;

  new_root_left->next  = new_root_right;
  new_root_right->next = new_root_left;

  pll_unode_t *left_child  = root_location.edge;
  pll_unode_t *right_child = left_child->back;

  left_child->back    = new_root_left;
  new_root_left->back = left_child;
  left_child->length = new_root_left->length = root_location.brlen();

  right_child->back    = new_root_right;
  new_root_right->back = right_child;
  right_child->length  = new_root_right->length =
      root_location.brlen_compliment();

  unsigned int total_unodes = _tree->inner_count * 3 + _tree->tip_count;

  _tree->inner_count += 1;
  _tree->edge_count += 1;
  _tree->vroot = new_root_left;

  new_root_left->clv_index = new_root_right->clv_index = tree_size - 1;
  new_root_left->scaler_index = new_root_right->scaler_index =
      static_cast<int>(_tree->inner_count - 1);

  new_root_left->pmatrix_index = left_child->pmatrix_index;

  new_root_right->node_index = total_unodes + 2;
  right_child->pmatrix_index = new_root_right->pmatrix_index =
      _tree->edge_count - 1;
  new_root_right->pmatrix_index = _tree->edge_count - 1;

  _current_rl = root_location;
  _rooted     = true;
}

void rooted_tree_t::update_root(root_location_t root) {
  if (root.edge != _tree->vroot) {
    throw std::runtime_error("Provided root doesn't match the current tree");
  }

  pll_unode_t *right_root = root.edge;
  pll_unode_t *left_root  = root.edge->next;

  right_root->length = right_root->back->length = root.brlen();
  left_root->length = left_root->back->length = root.brlen_compliment();
}

void rooted_tree_t::unroot() {
  pll_unode_t *left_child  = _tree->vroot->back;
  pll_unode_t *right_child = _tree->vroot->next->back;
  pll_unode_t *root_left   = _tree->vroot;
  pll_unode_t *root_right  = _tree->vroot->next;

  right_child->back = left_child;
  left_child->back  = right_child;

  right_child->length = left_child->length = _current_rl.saved_brlen;

  clean_root_unode(root_left);
  clean_root_unode(root_right);

  _tree->vroot = left_child->next != nullptr ? left_child : right_child;
  if (_tree->vroot->next == nullptr) {
    throw std::runtime_error("unrooted to a tip");
  }

  _tree->inner_count -= 1;
  _tree->edge_count -= 1;

  right_child->pmatrix_index = left_child->pmatrix_index;
  _rooted                    = false;
}

bool rooted_tree_t::rooted() const {
  return _tree->vroot->next->next == _tree->vroot;
}

std::tuple<std::vector<pll_operation_t>,
           std::vector<unsigned int>,
           std::vector<double>>
rooted_tree_t::generate_operations(const root_location_t &new_root) {
  root_by(new_root);
  auto trav_buf = full_traverse();
  debug_string(EMIT_LEVEL_DEBUG, "traversal after root");

  for (auto node : trav_buf) {
    debug_print(EMIT_LEVEL_DEBUG,
                "traversal node label: %s, pmatrix index: %d, clv index: %d",
                (node->label != nullptr ? node->label : "null"),
                node->pmatrix_index,
                node->clv_index);
  }

  std::vector<pll_operation_t> ops(trav_buf.size());
  std::vector<unsigned int>    pmatrix_indices(trav_buf.size());
  std::vector<double>          branch_lengths(trav_buf.size());

  unsigned int op_count     = 0;
  unsigned int matrix_count = 0;

  pll_utree_create_operations(trav_buf.data(),
                              static_cast<unsigned int>(trav_buf.size() - 1),
                              branch_lengths.data(),
                              pmatrix_indices.data(),
                              ops.data(),
                              &matrix_count,
                              &op_count);

  ops.resize(op_count + 1);
  pmatrix_indices.resize(matrix_count);
  branch_lengths.resize(matrix_count);

  auto         root_op_it         = ops.end() - 1;
  pll_unode_t *root_node          = *(trav_buf.end() - 1);
  root_op_it->parent_clv_index    = root_node->clv_index;
  root_op_it->parent_scaler_index = root_node->scaler_index;

  root_op_it->child1_clv_index    = root_node->back->clv_index;
  root_op_it->child1_scaler_index = root_node->back->scaler_index;
  root_op_it->child1_matrix_index = root_node->back->pmatrix_index;

  root_op_it->child2_clv_index    = root_node->next->back->clv_index;
  root_op_it->child2_scaler_index = root_node->next->back->scaler_index;
  root_op_it->child2_matrix_index = root_node->next->back->pmatrix_index;

  return std::make_tuple(ops, pmatrix_indices, branch_lengths);
}

std::tuple<pll_operation_t, std::vector<unsigned int>, std::vector<double>>
rooted_tree_t::generate_derivative_operations(const root_location_t &root) {
  root_by(root);

  pll_operation_t           op;
  std::vector<unsigned int> pmatrix_indices(2);
  std::vector<double>       branch_lengths(2);

  pll_unode_t *vroot = _tree->vroot;

  op.parent_clv_index    = root_clv_index();
  op.parent_scaler_index = root_scaler_index();

  op.child1_clv_index    = vroot->back->clv_index;
  op.child1_matrix_index = vroot->back->pmatrix_index;
  op.child1_scaler_index = vroot->back->scaler_index;
  pmatrix_indices[0]     = vroot->back->pmatrix_index;
  branch_lengths[0]      = vroot->back->length;

  op.child2_clv_index    = vroot->next->back->clv_index;
  op.child2_matrix_index = vroot->next->back->pmatrix_index;
  op.child2_scaler_index = vroot->next->back->scaler_index;
  pmatrix_indices[1]     = vroot->next->back->pmatrix_index;
  branch_lengths[1]      = vroot->next->back->length;

  return std::make_tuple(op, pmatrix_indices, branch_lengths);
}

std::string rooted_tree_t::newick(bool annotations) const {
  if (!_root_annotations.empty() && annotations) {
    auto fold_operation = [](std::string                         a,
                             std::pair<std::string, std::string> b) {
      return std::move(a) + ':' + std::move(b.first) + '='
             + std::move(b.second);
    };

    for (auto kv : _root_annotations) {
      if (kv.second.size() == 0) { continue; }
      debug_print(
          EMIT_LEVEL_DEBUG, "number of annotation: %lu", kv.second.size());
      auto        root_location = kv.first;
      std::string start         = "[&&NHX";
      std::string annotation    = std::accumulate(
          kv.second.begin(), kv.second.end(), start, fold_operation);
      annotation += ']';

      if (root_location->data) {
        debug_print(-1,
                    "root_location->data: %p, int: %ld",
                    root_location->data,
                    (int64_t)root_location->data);

        throw std::runtime_error("The node data was not empty when we tried to "
                                 "put node annotations in");
      }

      debug_print(
          EMIT_LEVEL_DEBUG, "calculated annotation: %s", annotation.c_str());
      char *new_label = (char *)calloc(sizeof(char), (annotation.size() + 1));
      memcpy(new_label, annotation.data(), sizeof(char) * annotation.size());
      debug_print(EMIT_LEVEL_DEBUG, "new label: %s", new_label);

      root_location->data = new_label;
    }
  }
  auto serialize_node = [](const pll_unode_t *n) {
    auto tmp = (n->label ? std::string(n->label) : std::string()) + ':'
               + std::to_string(n->length)
               + (n->data ? std::string((char *)n->data) : std::string());
    char *node_string = (char *)calloc(sizeof(char), tmp.size() + 1);
    memcpy(node_string, tmp.data(), tmp.size());
    return node_string;
  };
  char *newick_string = pll_utree_export_newick(_tree->vroot, serialize_node);
  std::string ret{newick_string};
  free(newick_string);
  return ret;
}

void rooted_tree_t::show_tree() const {
  pll_utree_show_ascii(_tree->vroot,
                       PLL_UTREE_SHOW_LABEL | PLL_UTREE_SHOW_BRANCH_LENGTH);
}

bool rooted_tree_t::branch_length_sanity_check() const {
  auto nodes = full_traverse();
  nodes.pop_back();
  std::sort(nodes.begin(), nodes.end(), [](pll_unode_t *a, pll_unode_t *b) {
    return a->length < b->length;
  });

  size_t median_index1 = (nodes.size() - 1) / 2;
  size_t median_index2 = (nodes.size()) / 2;
  double median =
      (nodes[median_index1]->length + nodes[median_index2]->length) / 2.0;

  if (median * 10.0 < (*(nodes.end() - 1))->length
      || ((*nodes.begin())->length) < median / 10.0) {
    return false;
  }
  return true;
}

bool rooted_tree_t::sanity_check() const {
  return branch_length_sanity_check();
}

static void tag_nodes(pll_unode_t *n) {
  pll_unode_t *start = n;
  do {
    n->data = (void *)0xdeadbeef;
    n       = n->next;
  } while (n != nullptr && n != start);
}

static void untag_nodes(pll_unode_t *n) {
  pll_unode_t *start = n;
  do {
    n->data = (void *)nullptr;
    n       = n->next;
  } while (n != nullptr && n != start);
}

void rooted_tree_t::find_path(pll_unode_t *n1, pll_unode_t *n2) {
  pll_unode_t *start   = n1;
  pll_unode_t *current = n1;
  do {
    if (find_path_recurse(current->back, n2)) break;
    current = current->next;
  } while (current != nullptr && current != start);
}

bool rooted_tree_t::find_path_recurse(pll_unode_t *n1, pll_unode_t *n2) {
  if (n1 == n2) {
    tag_nodes(n1);
    return true;
  }
  if (n1->next) {
    pll_unode_t *start = n1;
    n1                 = n1->next;
    while (start != n1) {
      if (n1 == n2) {
        tag_nodes(n1);
        return true;
      } else if (find_path_recurse(n1->back, n2)) {
        tag_nodes(n1);
        return true;
      }
      n1 = n1->next;
    }
  } else {
    if (n1 == n2) { tag_nodes(n1); }
    return n1 == n2;
  }
  return false;
}

std::tuple<std::vector<pll_operation_t>,
           std::vector<unsigned int>,
           std::vector<double>>
rooted_tree_t::generate_root_update_operations(
    const root_location_t &new_root) {

  if (new_root.edge == _current_rl.edge
      || new_root.edge == _current_rl.edge->back) {
    return {};
  }

  auto old_root = _current_rl;
  root_by(new_root);
  find_path(old_root.edge, _tree->vroot);
  /* cover the old root */
  tag_nodes(old_root.edge);
  tag_nodes(old_root.edge->back);

  tag_nodes(_tree->vroot->back);
  tag_nodes(_tree->vroot->next->back);
  std::vector<pll_unode_t *> trav_buf(inner_count() + tip_count());
  unsigned int               trav_size = 0;

  auto update_root_callback = [](pll_unode_t *n) -> int {
    if (n->data == (void *)0xdeadbeef) {
      n->data = nullptr;
      return PLL_SUCCESS;
    }
    return PLL_FAILURE;
  };

  pll_utree_traverse(_tree->vroot,
                     PLL_TREE_TRAVERSE_POSTORDER,
                     update_root_callback,
                     trav_buf.data(),
                     &trav_size);
  trav_buf.resize(trav_size);

  std::vector<pll_operation_t> ops(trav_buf.size());
  std::vector<unsigned int>    pmatrix_indices(trav_buf.size());
  std::vector<double>          branch_lengths(trav_buf.size());

  assert_string(trav_buf.size() != 0,
                "traversal buffer when updating the root had size zero");

  unsigned int op_count     = 0;
  unsigned int matrix_count = 0;

  pll_utree_create_operations(trav_buf.data(),
                              static_cast<unsigned int>(trav_buf.size() - 1),
                              branch_lengths.data(),
                              pmatrix_indices.data(),
                              ops.data(),
                              &matrix_count,
                              &op_count);

  ops.resize(op_count);
  pmatrix_indices.resize(matrix_count);
  branch_lengths.resize(matrix_count);

  ops.resize(op_count + 1);
  pmatrix_indices.resize(matrix_count);
  branch_lengths.resize(matrix_count);

  auto         root_op_it         = ops.end() - 1;
  pll_unode_t *root_node          = _tree->vroot;
  root_op_it->parent_clv_index    = root_node->clv_index;
  root_op_it->parent_scaler_index = root_node->scaler_index;

  root_op_it->child1_clv_index    = root_node->back->clv_index;
  root_op_it->child1_scaler_index = root_node->back->scaler_index;
  root_op_it->child1_matrix_index = root_node->back->pmatrix_index;

  root_op_it->child2_clv_index    = root_node->next->back->clv_index;
  root_op_it->child2_scaler_index = root_node->next->back->scaler_index;
  root_op_it->child2_matrix_index = root_node->next->back->pmatrix_index;

  clear_traversal_data();
  untag_nodes(old_root.edge);
  untag_nodes(old_root.edge->back);

  untag_nodes(_tree->vroot->back);
  untag_nodes(_tree->vroot->next->back);

  return std::make_tuple(ops, pmatrix_indices, branch_lengths);
}

void rooted_tree_t::clear_traversal_data() {
  for (unsigned int i = 0; i < _tree->tip_count; ++i) {
    _tree->nodes[i]->data = nullptr;
  }
  for (unsigned int i = _tree->tip_count;
       i < _tree->tip_count + _tree->inner_count;
       ++i) {
    pll_unode_t *start = _tree->nodes[i];
    pll_unode_t *node  = start;
    do {
      node->data = nullptr;
      node       = node->next;
    } while (node && node != start);
  }
  pll_unode_t *start = _tree->vroot;
  pll_unode_t *node  = start;
  do {
    node->data = nullptr;
    node       = node->next;
  } while (node && node != start);
}

root_location_t rooted_tree_t::current_root() const {
  if (!rooted())
    throw std::runtime_error("Failed to return root, tree is unrooted");
  return _current_rl;
}

const std::vector<root_location_t> &rooted_tree_t::roots() const {
  return _roots;
}

void rooted_tree_t::annotate_node(size_t             node_id,
                                  const std::string &key,
                                  const std::string &value) {
  annotate_node(_roots[node_id], key, value);
}

void rooted_tree_t::annotate_node(const root_location_t &node_index,
                                  const std::string &    key,
                                  const std::string &    value) {
  annotate_node(node_index.edge, key, value);
}

void rooted_tree_t::annotate_node(pll_unode_t *      node_id,
                                  const std::string &key,
                                  const std::string &value) {
  _root_annotations[node_id].emplace_back(key, value);
}

void rooted_tree_t::annotate_ratio(size_t node_id, double ratio) {
  annotate_ratio(_roots[node_id], ratio);
}

void rooted_tree_t::annotate_ratio(const root_location_t &node_index,
                                   double                 ratio) {
  annotate_branch(
      node_index, "alpha", std::to_string(ratio), std::to_string(1 - ratio));
}

void rooted_tree_t::annotate_lh(size_t node_index, double lh) {
  annotate_ratio(_roots[node_index], lh);
}

void rooted_tree_t::annotate_lh(const root_location_t &node_index, double lh) {
  annotate_branch(node_index, "LLH", std::to_string(lh));
}

void rooted_tree_t::annotate_branch(size_t             node_id,
                                    const std::string &key,
                                    const std::string &value) {
  annotate_branch(_roots[node_id], key, value);
}

void rooted_tree_t::annotate_branch(const root_location_t &rl,
                                    const std::string &    key,
                                    const std::string &    value) {
  annotate_branch(rl, key, value, value);
}

void rooted_tree_t::annotate_branch(const root_location_t &rl,
                                    const std::string &    key,
                                    const std::string &    left_value,
                                    const std::string &    right_value) {
  annotate_node(rl.edge, key, left_value);
  size_t       node_count = 0;
  pll_unode_t *start      = rl.edge->back;
  pll_unode_t *cur        = start;

  if (cur->next) {
    do {
      node_count++;
      cur = cur->next;
    } while (start != cur);
  } else {
    node_count = 1;
  }

  if (node_count > 2) {
    annotate_node(rl.edge->back, key, right_value);
  } else { // We have a root in this case (or a "fake node")
    annotate_node(rl.edge->back->next->back, key, right_value);
  }
}

std::unordered_map<pll_unode_t *, pll_unode_t *>
rooted_tree_t::make_node_bijection(const rooted_tree_t &other) {
  std::unordered_map<pll_unode_t *, pll_unode_t *> bijection;
  auto                                             our_nodes = full_traverse();
  auto their_nodes = other.full_traverse();

  if (our_nodes.size() != their_nodes.size()) {
    throw std::runtime_error{
        "Can't create a bijection between two different sized trees"};
  }

  for (size_t i = 0; i < our_nodes.size(); ++i) {
    pll_unode_t *our_cur_node   = our_nodes[i];
    pll_unode_t *their_cur_node = their_nodes[i];
    do {
      bijection[their_cur_node] = our_cur_node;
      our_cur_node              = our_cur_node->next;
      their_cur_node            = their_cur_node->next;
      if (!(our_cur_node == nullptr && their_cur_node == nullptr)
          && !(our_cur_node != nullptr && their_cur_node != nullptr)) {
        throw std::runtime_error{"There is a difference in topology discoverd "
                                 "when trying to make a bijective map"};
      }
    } while (our_cur_node != nullptr && their_cur_node != nullptr
             && our_cur_node != our_nodes[i]
             && their_cur_node != their_nodes[i]);
  }
  return bijection;
}

void rooted_tree_t::copy_annotations(const rooted_tree_t &other) {
  auto bijection = make_node_bijection(other);

  for (auto &kv : other._root_annotations) {
    auto &node                            = kv.first;
    auto &annotation                      = kv.second;
    _root_annotations[bijection.at(node)] = annotation;
  }
}
