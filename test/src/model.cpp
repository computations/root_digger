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
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    model_t model{tree, msa, {1}, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
  }
}

TEST_CASE("model_t constructor with partitions", "[model_t]") {
  auto &ds = data_files_dna["101.phy"];
  msa_t unparted_msa{ds.first};
  rooted_tree_t tree{ds.second};
  msa_partitions_t parts;
  parts.push_back(parse_partition_info("DNA, PART_0= 0-100"));
  parts.push_back(parse_partition_info("DNA, PART_1= 101-200"));
  auto msa = unparted_msa.partition(parts);
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1,1}, seed, false};
}

TEST_CASE("model_t compute lh", "[model_t]") {
  auto &ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) {
      continue;
    }
    double lh = model.compute_lh(tree.root_location(i));
    CHECK(std::isfinite(lh));
    CHECK(lh < 0.0);
  }
}

TEST_CASE("model_t compute lh, all data", "[!hide][all_data][model_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    model_t model{tree, msa, {1}, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
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

TEST_CASE("model_t compute dlh/da", "[model_t][opt]") {
  auto &ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
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

TEST_CASE("model_t compute dlh/da, all data",
          "[!hide][all_data][model_t][opt]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    model_t model{tree, msa, {1}, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
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

TEST_CASE("model_t optimize root locations on individual roots",
          "[model_t][opt]") {
  auto &ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) {
      continue;
    }
    model.compute_lh(tree.root_location(i));
    model.optimize_alpha(tree.root_location(i), 1e-7);
  }
}

TEST_CASE("model_t optimize root locations on individual roots, all data",
          "[!hide][all_data][model_t][opt]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    model_t model{tree, msa, {1}, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) {
        continue;
      }
      model.compute_lh(tree.root_location(i));
      model.optimize_alpha(tree.root_location(i), 1e-7);
    }
  }
}

TEST_CASE("model_t optimize root locations with beg points", "[model_t][opt]") {
  auto &ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) {
      continue;
    }
    auto rl = tree.root_location(i);
    rl.brlen_ratio = 0.0;
    model.compute_lh(rl);
    model.optimize_alpha(rl, 1e-7);
  }
}

TEST_CASE("model_t optimize root locations with beg points, all data",
          "[!hide][all_data][model_t][opt]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    model_t model{tree, msa, {1}, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) {
        continue;
      }
      auto rl = tree.root_location(i);
      rl.brlen_ratio = 0.0;
      model.compute_lh(rl);
      model.optimize_alpha(rl, 1e-7);
    }
  }
}

TEST_CASE("model_t optimize root locations with end points", "[model_t][opt]") {
  auto &ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) {
      continue;
    }
    auto rl = tree.root_location(i);
    rl.brlen_ratio = 1.0;
    model.compute_lh(rl);
    model.optimize_alpha(rl, 1e-7);
  }
}

TEST_CASE("model_t optimize root locations with end points, all data",
          "[!hide][all_data][model_t][opt]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    model_t model{tree, msa, {1}, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) {
        continue;
      }
      auto rl = tree.root_location(i);
      rl.brlen_ratio = 1.0;
      model.compute_lh(rl);
      model.optimize_alpha(rl, 1e-7);
    }
  }
}

TEST_CASE("model_t optimize whole tree", "[model_t]") {
  auto &ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  auto rl = model.optimize_root_location().first;
  CHECK(rl.brlen_ratio >= 0.0);
  CHECK(rl.brlen_ratio <= 1.0);
}

TEST_CASE("model_t optimize whole tree, all data",
          "[!hide][all_data][model_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    model_t model{tree, msa, {1}, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    auto rl = model.optimize_root_location().first;
    CHECK(rl.brlen_ratio >= 0.0);
    CHECK(rl.brlen_ratio <= 1.0);
  }
}

TEST_CASE("model_t liklihood computation", "[model_t]") {
  auto &ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
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

TEST_CASE("model_t liklihood computation, all data",
          "[!hide][all_data][model_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t seed = std::rand();
    model_t model{tree, msa, {1}, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
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

TEST_CASE("model_t optimize all", "[model_t]") {
  auto ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  auto initial_rl = model.optimize_root_location();
  auto tmp = model.optimize_all(1, 0.0, 1e-3, 1e-3, 1e-3, 1e12);
  auto final_rl = tmp.first;
  auto final_lh = tmp.second;
  CHECK(final_lh >= initial_rl.second);
  CHECK(model.compute_lh(final_rl) == Approx(final_lh));
}

TEST_CASE("model_t optimize all, slow", "[!hide][all_data][model_t]") {
  auto ds = data_files_dna["101.phy"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  auto initial_rl = model.optimize_root_location();
  auto tmp = model.optimize_all(1, 0.0, 1e-7, 1e-7, 1e-7, 1e7);
  auto final_rl = tmp.first;
  auto final_lh = tmp.second;
  CHECK(final_lh >= initial_rl.second);
  CHECK(model.compute_lh(final_rl) == Approx(final_lh));
}

TEST_CASE("model_t move root", "[!hide][model_t]") {
  auto ds = data_files_dna["101.phy"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, {1}, seed, false};

  model.set_subst_rates(
      0, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
  model.set_freqs(0, {0.25, 0.25, 0.25, 0.25});

  auto root_lh = model.compute_all_root_lh();
  for (size_t i = 0; i < root_lh.size(); ++i) {
    for (size_t j = i + 1; j < root_lh.size(); ++j) {
      CHECK(root_lh[i] == Approx(root_lh[j]));
    }
  }
}
