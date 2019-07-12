#ifndef __RD_TREE_HPP_
#define __RD_TREE_HPP_

extern "C" {
#include <libpll/pll.h>
}
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#define GENERATE_AND_UNPACK_OPS(TREE, RL, OPS, PM, BR)                         \
  {                                                                            \
    auto results = TREE.generate_operations(RL);                               \
    OPS = std::move(std::get<0>(results));                                     \
    PM = std::move(std::get<1>(results));                                      \
    BR = std::move(std::get<2>(results));                                      \
  }

struct root_location_t {
  pll_unode_t *edge;
  double saved_brlen;
  double brlen_ratio;
  constexpr inline double brlen() const { return saved_brlen * brlen_ratio; }
  constexpr inline double brlen_compliment() const {
    return saved_brlen * (1 - brlen_ratio);
  }
  std::string label() const{
    return edge->label != nullptr ? edge->label : "(null)";
  }
};

pll_utree_t *parse_tree_file(const std::string &tree_filename);

class rooted_tree_t {
public:
  rooted_tree_t() : _tree{nullptr}, _roots{} {}

  rooted_tree_t(const std::string &tree_filename)
      : _tree{parse_tree_file(tree_filename)}, _roots{} {
    if (_tree == nullptr) {
      throw std::invalid_argument("Tree file could not be parsed");
    }
    generate_root_locations();
  }

  rooted_tree_t(rooted_tree_t &&other)
      : _tree{std::move(other._tree)}, _roots{std::move(other._roots)} {
    other._tree = nullptr;
  }

  rooted_tree_t(const rooted_tree_t &);
  ~rooted_tree_t();

  rooted_tree_t &operator=(rooted_tree_t &&);
  rooted_tree_t &operator=(const rooted_tree_t &);

  root_location_t root_location(size_t) const;
  size_t root_count() const;

  unsigned int tip_count() const;
  unsigned int inner_count() const;
  unsigned int branch_count() const;

  unsigned int root_clv_index() const;
  unsigned int root_scaler_index() const;

  std::unordered_map<std::string, unsigned int> label_map() const;
  std::unordered_set<std::string> label_set() const;

  std::tuple<std::vector<pll_operation_t>, std::vector<unsigned int>,
             std::vector<double>>
  generate_operations(const root_location_t &);

  std::tuple<pll_operation_t, std::vector<unsigned int>, std::vector<double>>
  generate_derivative_operations(const root_location_t &root);

  void root_by(const root_location_t &);
  void update_root(root_location_t);
  void unroot();
  bool rooted() const;
  bool branch_length_sanity_check() const;
  bool sanity_check() const;

  std::string newick() const;

  void show_tree() const;

private:
  void generate_root_locations();
  void deduplicate_roots();
  std::vector<pll_unode_t *> full_traverse() const;
  std::vector<pll_unode_t *> edge_traverse() const;

  pll_utree_t *_tree;
  std::vector<root_location_t> _roots;
};

#endif
