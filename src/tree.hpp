#ifndef __RD_TREE_HPP_
#define __RD_TREE_HPP_

extern "C" {
#include <libpll/pll.h>
}
#include <string>

pll_utree_t *parse_tree_file(const std::string &tree_filename);

#endif
