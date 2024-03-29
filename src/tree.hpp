#ifndef RD_TREE_HPP_
#define RD_TREE_HPP_

extern "C" {
#include <corax/corax.h>
}
#include "debug.h"
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#define GENERATE_AND_UNPACK_OPS(TREE, RL, OPS, PM, BR)                         \
  {                                                                            \
    auto results = TREE.generate_operations(RL);                               \
    OPS          = std::move(std::get<0>(results));                            \
    PM           = std::move(std::get<1>(results));                            \
    BR           = std::move(std::get<2>(results));                            \
  }

struct root_location_t {
  corax_unode_t *edge;
  size_t         id;
  double         saved_brlen;
  double         brlen_ratio;

  constexpr inline double brlen() const { return saved_brlen * brlen_ratio; }
  constexpr inline double brlen_compliment() const {
    return saved_brlen * (1 - brlen_ratio);
  }

  std::string label() const {
    return edge->label != nullptr ? edge->label : "(null)";
  }

  bool is_internal() const {
    return edge->next != nullptr && edge->back->next != nullptr;
  }
  bool is_external() const { return !is_internal(); }

  bool operator==(const root_location_t &other) const {
    return edge == other.edge && brlen_ratio == other.brlen_ratio;
  }
  bool operator!=(const root_location_t &other) const {
    return !(*this == other);
  }
};

corax_utree_t *parse_tree_file(const std::string &tree_filename);

class rooted_tree_t {
public:
  rooted_tree_t() : _tree{nullptr}, _roots{}, _rooted{false} {}

  rooted_tree_t(const std::string &tree_filename) :
      _tree{parse_tree_file(tree_filename)}, _roots{}, _rooted{false} {
    if (_tree == nullptr) {
      throw std::invalid_argument("Tree file could not be parsed");
    }
    generate_root_locations();
    add_root_space();
    sort_root_locations();
  }

  rooted_tree_t(rooted_tree_t &&other) :
      _tree{std::move(other._tree)},
      _roots{std::move(other._roots)},
      _root_annotations{std::move(other._root_annotations)},
      _rooted{other._rooted} {
    other._tree = nullptr;
  }

  rooted_tree_t(const rooted_tree_t &);
  ~rooted_tree_t();

  rooted_tree_t &operator=(rooted_tree_t &&);
  rooted_tree_t &operator=(const rooted_tree_t &);

  root_location_t root_location(size_t) const;
  root_location_t root_location(const std::string &) const;
  root_location_t root_location() const;

  root_location_t midpoint() const;

  std::vector<root_location_t> rank_midpoints() const;
  std::vector<root_location_t> rank_modified_mad() const;

  size_t       root_count() const;
  unsigned int tip_count() const;
  unsigned int inner_count() const;
  unsigned int branch_count() const;

  unsigned int root_clv_index() const;
  int          root_scaler_index() const;

  root_location_t                     current_root() const;
  const std::vector<root_location_t> &roots() const;

  std::vector<root_location_t> internal_root_locations() const;
  std::vector<root_location_t> external_root_locations() const;

  std::unordered_map<std::string, unsigned int> label_map() const;
  std::unordered_set<std::string>               label_set() const;

  std::tuple<std::vector<corax_operation_t>,
             std::vector<unsigned int>,
             std::vector<double>>
  generate_operations(const root_location_t &);

  std::tuple<corax_operation_t, std::vector<unsigned int>, std::vector<double>>
  generate_derivative_operations(const root_location_t &root);

  std::tuple<std::vector<corax_operation_t>,
             std::vector<unsigned int>,
             std::vector<double>>
  generate_root_update_operations(const root_location_t &new_root);

  void root_by(unsigned int root_id);
  void root_by(const root_location_t &);
  void update_root(root_location_t);
  void unroot();
  bool rooted() const;
  bool branch_length_sanity_check() const;
  bool sanity_check() const;

  std::string newick(bool annotations = true) const;

  void clear_newick_annotations();

  void show_tree() const;

  void annotate_node(const root_location_t &rl,
                     const std::string     &key,
                     const std::string     &value);
  void annotate_node(size_t             node_id,
                     const std::string &key,
                     const std::string &value);

  void annotate_branch(size_t             node_id,
                       const std::string &key,
                       const std::string &value);
  void annotate_branch(const root_location_t &rl,
                       const std::string     &key,
                       const std::string     &value);
  void annotate_branch(const root_location_t &rl,
                       const std::string     &key,
                       const std::string     &left_value,
                       const std::string     &right_value);

  void annotate_lh(size_t node_index, double lh);
  void annotate_lh(const root_location_t &node_index, double lh);

  void annotate_ratio(size_t node_id, double ratio);
  void annotate_ratio(const root_location_t &node_index, double ratio);

  std::vector<double> get_children_distance(corax_unode_t *rl);

  std::vector<std::pair<root_location_t, double>>
  apply_foreach_branch_map_reduce(
      const std::function<double(double, double, double)>      &map_func,
      const std::function<double(const std::vector<double> &)> &reduce_func)
      const;

private:
  void sort_root_locations();
  void generate_root_locations();
  void copy_root_locations(const rooted_tree_t &);
  void add_root_space();

  std::vector<corax_unode_t *> full_traverse() const;
  std::vector<corax_unode_t *> edge_traverse() const;

  void find_path(corax_unode_t *n1, corax_unode_t *n2);
  bool find_path_recurse(corax_unode_t *n1, corax_unode_t *n2);

  void clear_traversal_data();
  void clear_traversal_data(corax_unode_t *);

  void annotate_node(corax_unode_t     *node_id,
                     const std::string &key,
                     const std::string &value);

  void copy_annotations(const rooted_tree_t &other);

  std::unordered_map<corax_unode_t *, corax_unode_t *>
  make_node_bijection(const rooted_tree_t &other);

  std::vector<double> get_forward_children_distance(corax_unode_t *rl) const;
  std::vector<double> get_backward_children_distance(corax_unode_t *rl) const;

  corax_utree_t               *_tree;
  root_location_t              _current_rl;
  std::vector<root_location_t> _roots;
  std::unordered_map<corax_unode_t *,
                     std::vector<std::pair<std::string, std::string>>>
       _root_annotations;
  bool _rooted;
};

#endif
