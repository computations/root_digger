#ifndef __RD_TREE_HPP_
#define __RD_TREE_HPP_

extern "C" {
#include <libpll/pll.h>
}
#include <stdexcept>
#include <string>
#include <unordered_map>
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
  double brlen_ratio;
  constexpr inline double brlen() const { return edge->length * brlen_ratio; }
  constexpr inline double brlen_compliment() const {
    return edge->length * (1 - brlen_ratio);
  }
};

pll_utree_t *parse_tree_file(const std::string &tree_filename);

class rooted_tree_t {
public:
  rooted_tree_t()
      : _tree{nullptr}, _roots{}, _root_clv_index{0}, _root_scaler_index{0} {}
  rooted_tree_t(const std::string &tree_filename)
      : _tree{parse_tree_file(tree_filename)}, _roots{},
        _root_clv_index{0},
        _root_scaler_index{0} {
    if (_tree == nullptr) {
      throw std::invalid_argument("Tree file could not be parsed");
    }
    _root_clv_index = _tree->inner_count+ _tree->tip_count;
    _root_scaler_index = _tree->inner_count+ _tree->tip_count;
    generate_root_locations();
  }
  rooted_tree_t(rooted_tree_t &&other)
      : _tree{std::move(other._tree)}, _roots{std::move(other._roots)},
        _root_clv_index{other._root_clv_index}, _root_scaler_index{
                                                    other._root_scaler_index} {
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
  unsigned int branches() const;

  unsigned int root_clv_index() const;
  unsigned int root_scaler_index() const;

  std::unordered_map<std::string, unsigned int> label_map() const;

  std::tuple<std::vector<pll_operation_t>, std::vector<unsigned int>,
             std::vector<double>>
  generate_operations(const root_location_t &) const;

private:
  std::vector<pll_unode_t *> full_tree_traverse() const;
  std::vector<pll_unode_t *> full_tree_traverse(pll_unode_t *) const;
  void generate_root_locations();

  pll_utree_t *_tree;
  std::vector<root_location_t> _roots;
  unsigned int _root_clv_index;
  unsigned int _root_scaler_index;
};

#endif
