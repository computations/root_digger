#include "data.hpp"
#include <catch2/catch.hpp>
#include <cmath>
#include <debug.h>
#include <model.hpp>

model_params_t params[] = {
    {1, 2.5, 1, 1, 1, 2.5, 2.5, 1, 1, 1, 2.5, 1},
    {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
    {1.0, 1.2, 1.3, 1.2, 1.0, 1.3, 1.2, 1.1, 1.1, 1.4, 1.0, 1.0},
    {.34, .42, .24, .74, .16, .88, .75, .54, .20, .06, .08, .41},
};

model_params_t freqs[] = {
    {.25, .25, .25, .25},
};

TEST_CASE("model_t constructor", "[model_t]") {
  for (auto &mp : params) {
    for (auto &f : freqs) {
      for (auto &ds : data_files_dna) {
        msa_t msa{ds.first};
        rooted_tree_t tree{ds.second};
        uint64_t seed = std::rand();
        model_t{mp, tree, msa, f, seed};
      }
    }
  }
}

TEST_CASE("model_t constructor, emperical freqs", "[model_t]") {
  for (auto &mp : params) {
    for (auto &ds : data_files_dna) {
      msa_t msa{ds.first};
      rooted_tree_t tree{ds.second};
      uint64_t seed = std::rand();
      if (msa.length() == 1) {
        CHECK_THROWS(model_t{mp, tree, msa, seed});
      } else {
        CHECK_NOTHROW(model_t{mp, tree, msa, seed});
      }
    }
  }
}

TEST_CASE("model_t constructor, random rates", "[model_t]") {
  for (auto &ds : data_files_dna) {
    msa_t msa{ds.first};
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    if (msa.length() == 1) {
      CHECK_THROWS(model_t{tree, msa, seed});
    } else {
      CHECK_NOTHROW(model_t{tree, msa, seed});
    }
  }
}

TEST_CASE("model_t compute lh", "[model_t]") {
  for (auto &ds : data_files_dna) {
    for (auto &f : freqs) {
      for (auto &mp : params) {
        msa_t msa{ds.first};
        rooted_tree_t tree{ds.second};
        uint64_t seed = std::rand();
        model_t model{mp, tree, msa, f, seed};
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
}

TEST_CASE("model_t compute dlh/da", "[model_t][opt]") {
  for (auto &ds : data_files_dna) {
    for (auto &f : freqs) {
      for (auto &mp : params) {
        msa_t msa{ds.first};
        rooted_tree_t tree{ds.second};
        uint64_t seed = std::rand();
        model_t model{mp, tree, msa, f, seed};
        for (size_t i = 0; i < tree.root_count(); ++i) {
          if (tree.root_location(i).edge->back == nullptr) {
            continue;
          }
          model.compute_lh(tree.root_location(i));
          auto dlh = model.compute_dlh(tree.root_location(i));
          CHECK(std::isfinite(dlh.lh));
          CHECK(std::isfinite(dlh.dlh));
        }
      }
    }
  }
}

TEST_CASE("model_t optimize root locations on individual roots",
          "[model_t][opt]") {
  for (auto &ds : data_files_dna) {
    for (auto &f : freqs) {
      for (auto &mp : params) {
        msa_t msa{ds.first};
        rooted_tree_t tree{ds.second};
        uint64_t seed = std::rand();
        model_t model{mp, tree, msa, f, seed};
        for (size_t i = 0; i < tree.root_count(); ++i) {
          if (tree.root_location(i).edge->back == nullptr) {
            continue;
          }
          model.optimize_alpha(tree.root_location(i));
        }
      }
    }
  }
}

TEST_CASE("model_t optimize root locations with beg points", "[model_t][opt]") {
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      for (auto &f : freqs) {
        msa_t msa{ds.first};
        rooted_tree_t tree{ds.second};
        uint64_t seed = std::rand();
        model_t model{mp, tree, msa, f, seed};
        for (size_t i = 0; i < tree.root_count(); ++i) {
          if (tree.root_location(i).edge->back == nullptr) {
            continue;
          }
          auto rl = tree.root_location(i);
          rl.brlen_ratio = 0.0;
          model.optimize_alpha(rl);
        }
      }
    }
  }
}

TEST_CASE("model_t optimize root locations with end points", "[model_t][opt]") {
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      for (auto &f : freqs) {
        msa_t msa{ds.first};
        rooted_tree_t tree{ds.second};
        uint64_t seed = std::rand();
        model_t model{mp, tree, msa, f, seed};
        for (size_t i = 0; i < tree.root_count(); ++i) {
          if (tree.root_location(i).edge->back == nullptr) {
            continue;
          }
          auto rl = tree.root_location(i);
          rl.brlen_ratio = 1.0;
          model.optimize_alpha(rl);
        }
      }
    }
  }
}

TEST_CASE("model_t optimize whole tree", "[model_t]") {
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      for (auto &f : freqs) {
        msa_t msa{ds.first};
        rooted_tree_t tree{ds.second};
        uint64_t seed = std::rand();
        model_t model{mp, tree, msa, f, seed};
        auto rl = model.optimize_root_location().first;
        CHECK(rl.brlen_ratio >= 0.0);
        CHECK(rl.brlen_ratio <= 1.0);
      }
    }
  }
}

TEST_CASE("model_t liklihood computation", "[model_t]") {
  for (auto &ds : data_files_dna) {
    for (auto &mp : params) {
      for (auto &f : freqs) {
        msa_t msa{ds.first};
        rooted_tree_t tree{ds.second};
        uint64_t seed = std::rand();
        model_t model{mp, tree, msa, f, seed};
        for (size_t i = 0; i < tree.root_count(); ++i) {
          if (tree.root_location(i).edge->back == nullptr) {
            continue;
          }
          auto rl = tree.root_location(i);
          model.compute_lh(rl);
          CHECK(std::fabs(model.compute_lh(rl) - model.compute_lh_root(rl)) ==
                Approx(0.0));
        }
      }
    }
  }
}

TEST_CASE("model_t optimize all", "[model_t]"){
  auto ds = data_files_dna[1];
  msa_t msa{ds.first};
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa,seed};
  auto initial_rl = model.optimize_root_location();
  auto final_rl = model.optimize_all();
  CHECK(model.calculate_lh(final_rl) >= initial_rl.second);
}
