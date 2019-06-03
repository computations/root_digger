#include "data.hpp"
#include <catch2/catch.hpp>
#include <cmath>
#include <model.hpp>

model_params_t params[] = {
    {1, 2.5, 1, 1, 1, 2.5, 2.5, 1, 1, 1, 2.5, 1},
    {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
    {1.0, 1.2, 1.3, 1.2, 1.0, 1.3, 1.2, 1.1, 1.1, 1.4, 1.0, 1.0},
    {.34, .42, .24, .74, .16, .88, .75, .54, .20, .06, .08, .41},
};

TEST_CASE("model_t constructor", "[model_t]") {
  for (auto &mp : params) {
    for (auto &ds : data_files_dna) {
      msa_t msa{ds.first};
      rooted_tree_t tree{ds.second};
      model_t{mp, tree, msa};
    }
  }
}

TEST_CASE("model_t compute lh", "[model_t]") {
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      msa_t msa{ds.first};
      rooted_tree_t tree{ds.second};
      model_t model{mp, tree, msa};
      for (size_t i = 0; i < tree.root_count(); ++i) {
        if (tree.root_location(i).edge->back == nullptr) {
          continue;
        }
        double lh = model.compute_lh(tree.root_location(i));
        CHECK(std::isfinite(lh));
        CHECK(lh < 0.0);
      }
    }
  }
}

TEST_CASE("model_t compute dlh/da", "[model_t][opt]") {
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      msa_t msa{ds.first};
      rooted_tree_t tree{ds.second};
      model_t model{mp, tree, msa};
      for (size_t i = 0; i < tree.root_count(); ++i) {
        if (tree.root_location(i).edge->back == nullptr) {
          continue;
        }
        double dlh = model.compute_dlh(tree.root_location(i));
        CHECK(std::isfinite(dlh));
      }
    }
  }
}

TEST_CASE("model_t optimize root locations on individual roots",
          "[model_t][opt]") {
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      msa_t msa{ds.first};
      rooted_tree_t tree{ds.second};
      model_t model{mp, tree, msa};
      for (size_t i = 0; i < tree.root_count(); ++i) {
        if (tree.root_location(i).edge->back == nullptr) {
          continue;
        }
        model.optimize_alpha(tree.root_location(i));
      }
    }
  }
}

TEST_CASE("model_t optimize root locations with extreme points",
          "[model_t][opt]") {
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      msa_t msa{ds.first};
      rooted_tree_t tree{ds.second};
      model_t model{mp, tree, msa};
      for (size_t i = 0; i < tree.root_count(); ++i) {
        if (tree.root_location(i).edge->back == nullptr) {
          continue;
        }
        auto rl = tree.root_location(i);
        rl.brlen_ratio = 0.0;
        model.optimize_alpha(tree.root_location(i));
      }
    }
  }
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      msa_t msa{ds.first};
      rooted_tree_t tree{ds.second};
      model_t model{mp, tree, msa};
      for (size_t i = 0; i < tree.root_count(); ++i) {
        if (tree.root_location(i).edge->back == nullptr) {
          continue;
        }
        auto rl = tree.root_location(i);
        rl.brlen_ratio = 1.0;
        model.optimize_alpha(tree.root_location(i));
      }
    }
  }
}
