#include "data.hpp"
#include <catch2/catch.hpp>
#include <cmath>
#include <model.hpp>

TEST_CASE("model_t constructor", "[model_t]") {
  model_params_t params{1.0, 1.2, 1.3, 1.2, 1.0, 1.3,
                        1.2, 1.1, 1.1, 1.4, 1.0, 1.0};
  for (auto &ds : data_files_dna) {
    msa_t msa{ds.first};
    rooted_tree_t tree{ds.second};
    model_t{params, tree, msa};
  }
}

/*
TEST_CASE("model_t compute lhe", "[model_t]") {
  model_params_t params{.34, .42, .24, .74, .16, .88,
                        .75, .54, .20, .06, .08, .41};
  for (auto &ds : data_files_dna) {
    msa_t msa{ds.first};
    rooted_tree_t tree{ds.second};
    model_t model{params, tree, msa};
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) {
        continue;
      }
      double lh = model.compute_lh(tree.root_location(i));
      CHECK(std::isfinite(lh));
      CHECK(lh < 0.0);
      std::cout << std::setprecision(100) << lh << std::endl;
    }
  }
}
*/
