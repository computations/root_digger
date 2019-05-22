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

struct root_location_t {
  pll_unode_t *edge;
  double brlen_ratio;
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
  }
  rooted_tree_t(rooted_tree_t &&other)
      : _tree{std::move(other._tree)}, _roots{std::move(other._roots)} {
    other._tree = nullptr;
  }
  rooted_tree_t(const rooted_tree_t &);
  ~rooted_tree_t();

  rooted_tree_t &operator=(rooted_tree_t &&);
  rooted_tree_t &operator=(const rooted_tree_t &);

  unsigned int tip_count() const;
  unsigned int inner_count() const;
  unsigned int branches() const;

  std::unordered_map<std::string, unsigned int> label_map() const;

  std::vector<pll_operation_t>
  generate_operations(const root_location_t &) const;

private:
  std::vector<pll_unode_t*> full_tree_traverse() const;
  std::vector<pll_unode_t*> full_tree_traverse(pll_unode_t*) const;
  void generate_root_locations();

  pll_utree_t *_tree;
  std::vector<root_location_t> _roots;
};

#endif
