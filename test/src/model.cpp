#include "data.hpp"
#include <catch2/catch.hpp>
#include <model.hpp>

TEST_CASE("model_t constructor", "[model_t]") {
  model_params_t params{1.0, 1.2, 1.3, 1.2, 1.0, 1.3,
                        1.2, 1.1, 1.1, 1.4, 1.0, 1.0};
  for (auto &ds : data_files_dna) {
    msa_t msa {ds.first};
    rooted_tree_t tree{ds.second};
    model_t {params, tree, msa};
  }
}
