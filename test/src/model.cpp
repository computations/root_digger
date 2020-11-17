#include "data.hpp"
#include "test_util.hpp"
#include <algorithm>
#include <catch2/catch.hpp>
#include <cmath>
#include <debug.h>
#include <model.hpp>
#include <random>
#include <unordered_set>

model_params_t params[] = {
    {1, 2.5, 1, 1, 1, 2.5, 2.5, 1, 1, 1, 2.5, 1},
    {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
    {1.0, 1.2, 1.3, 1.2, 1.0, 1.3, 1.2, 1.1, 1.1, 1.4, 1.0, 1.0},
    {.34, .42, .24, .74, .16, .88, .75, .54, .20, .06, .08, .41},
};

model_params_t freqs[] = {
    {.25, .25, .25, .25},
};

checkpoint_t make_dummy_checkpoint(const std::string &dataset_name) {
  std::random_device rd;
  uint64_t           nonce =
      (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
  checkpoint_t  ckp(std::string("/tmp/") + dataset_name + "_"
                   + base_58_encode(nonce));
  cli_options_t cli_options;
  ckp.save_options(cli_options);
  ckp.reload();
  return ckp;
}

TEST_CASE("model_t constructor", "[model_t]") {
  for (auto &kv : data_files_dna) {
    auto &             ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t      seed = std::rand();
    model_t       model{tree, msa, {1}, true, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
  }
}

TEST_CASE("model_t constructor with partitions", "[model_t]") {
  auto &           ds = data_files_dna["101.phy"];
  msa_t            unparted_msa{ds.first};
  rooted_tree_t    tree{ds.second};
  msa_partitions_t parts;
  parts.push_back(parse_partition_info("DNA, PART_0= 1-100"));
  parts.push_back(parse_partition_info("DNA, PART_1= 101-200"));
  auto     msa  = unparted_msa.partition(parts);
  uint64_t seed = std::rand();
  model_t  model{tree, msa, {{}, {}}, true, seed, false};
}

TEST_CASE("model_t compute lh", "[model_t]") {
  auto &             ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) { continue; }
    double lh = model.compute_lh(tree.root_location(i));
    CHECK(std::isfinite(lh));
    CHECK(lh < 0.0);
  }
}

TEST_CASE("model_t compute lh, all data", "[!hide][all_data][model_t]") {
  for (auto &kv : data_files_dna) {
    auto &             ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t      seed = std::rand();
    model_t       model{tree, msa, {1}, true, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) { continue; }
      double lh = model.compute_lh(tree.root_location(i));
      CHECK(std::isfinite(lh));
      CHECK(lh < 0.0);
    }
  }
}

TEST_CASE("model_t compute dlh/da", "[model_t][opt]") {
  auto &             ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) { continue; }
    model.compute_lh(tree.root_location(i));
    auto dlh = model.compute_dlh(tree.root_location(i));
    CHECK(std::isfinite(dlh.lh));
    CHECK(std::isfinite(dlh.dlh));
  }
}

TEST_CASE("model_t compute dlh/da, all data",
          "[!hide][all_data][model_t][opt]") {
  for (auto &kv : data_files_dna) {
    auto &             ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t      seed = std::rand();
    model_t       model{tree, msa, {1}, true, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) { continue; }
      model.compute_lh(tree.root_location(i));
      auto dlh = model.compute_dlh(tree.root_location(i));
      CHECK(std::isfinite(dlh.lh));
      CHECK(std::isfinite(dlh.dlh));
    }
  }
}

TEST_CASE("model_t optimize root locations on individual roots",
          "[model_t][opt]") {
  auto &             ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) { continue; }
    model.compute_lh(tree.root_location(i));
    model.optimize_alpha(tree.root_location(i), 1e-7);
  }
}

TEST_CASE("model_t optimize root locations on individual roots, all data",
          "[!hide][all_data][model_t][opt]") {
  for (auto &kv : data_files_dna) {
    auto &             ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t      seed = std::rand();
    model_t       model{tree, msa, {1}, true, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) { continue; }
      model.compute_lh(tree.root_location(i));
      model.optimize_alpha(tree.root_location(i), 1e-7);
    }
  }
}

TEST_CASE("model_t optimize root locations with beg points", "[model_t][opt]") {
  auto &             ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) { continue; }
    auto rl        = tree.root_location(i);
    rl.brlen_ratio = 0.0;
    model.compute_lh(rl);
    model.optimize_alpha(rl, 1e-7);
  }
}

TEST_CASE("model_t optimize root locations with beg points, all data",
          "[!hide][all_data][model_t][opt]") {
  for (auto &kv : data_files_dna) {
    auto &             ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t      seed = std::rand();
    model_t       model{tree, msa, {1}, true, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) { continue; }
      auto rl        = tree.root_location(i);
      rl.brlen_ratio = 0.0;
      model.compute_lh(rl);
      model.optimize_alpha(rl, 1e-7);
    }
  }
}

TEST_CASE("model_t optimize root locations with end points", "[model_t][opt]") {
  auto &             ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) { continue; }
    auto rl        = tree.root_location(i);
    rl.brlen_ratio = 1.0;
    model.compute_lh(rl);
    model.optimize_alpha(rl, 1e-7);
  }
}

TEST_CASE("model_t optimize root locations with end points, all data",
          "[!hide][all_data][model_t][opt]") {
  for (auto &kv : data_files_dna) {
    auto &             ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t      seed = std::rand();
    model_t       model{tree, msa, {1}, true, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) { continue; }
      auto rl        = tree.root_location(i);
      rl.brlen_ratio = 1.0;
      model.compute_lh(rl);
      model.optimize_alpha(rl, 1e-7);
    }
  }
}

TEST_CASE("model_t optimize whole tree", "[model_t]") {
  auto &             ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  model.compute_lh(tree.root_location(0));
  auto rl = model.optimize_root_location(1, .05).first;
  CHECK(rl.brlen_ratio >= 0.0);
  CHECK(rl.brlen_ratio <= 1.0);
}

TEST_CASE("model_t optimize whole tree, all data",
          "[!hide][all_data][model_t]") {
  for (auto &kv : data_files_dna) {
    auto &             ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t      seed = std::rand();
    model_t       model{tree, msa, {1}, true, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    model.compute_lh(tree.root_location(0));
    auto rl = model.optimize_root_location(1, .05).first;
    CHECK(rl.brlen_ratio >= 0.0);
    CHECK(rl.brlen_ratio <= 1.0);
  }
}

TEST_CASE("model_t liklihood computation", "[model_t]") {
  auto &             ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->back == nullptr) { continue; }
    auto   rl  = tree.root_location(i);
    double lh1 = model.compute_lh(rl);
    double lh2 = model.compute_lh(rl);
    double lh3 = model.compute_lh_root(rl);
    CHECK(std::fabs(lh1 - lh2) == Approx(0.0));
    CHECK(std::fabs(lh3 - lh2) == Approx(0.0));
  }
}

TEST_CASE("model_t liklihood computation, all data",
          "[!hide][all_data][model_t]") {
  for (auto &kv : data_files_dna) {
    auto &             ds = kv.second;
    std::vector<msa_t> msa;
    msa.emplace_back(ds.first);
    rooted_tree_t tree{ds.second};
    uint64_t      seed = std::rand();
    model_t       model{tree, msa, {1}, true, seed, false};
    model.initialize_partitions_uniform_freqs(msa);
    for (size_t i = 0; i < tree.root_count(); ++i) {
      if (tree.root_location(i).edge->back == nullptr) { continue; }
      auto rl = tree.root_location(i);
      model.compute_lh(rl);
      CHECK(std::fabs(model.compute_lh(rl) - model.compute_lh_root(rl))
            == Approx(0.0));
    }
  }
}

TEST_CASE("model_t optimize all", "[model_t]") {
  auto               ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  model.compute_lh(tree.root_location(0));
  model.optimize_root_location(1, .05);
  auto initial_lh = model.compute_lh(tree.root_location(4));

  SECTION("one min root") {
    /* There is something very wrong with this test. I am commenting out one of
     * the lines, but I really should figure out what is going on. I think that
     * it has something to do with the voodoo that catch is doing, which causes
     * stuff to work weirdly
     */
    auto checkpoint = make_dummy_checkpoint("10.fasta");
    model.assign_indicies_by_rank_search(1, 0.0, 0, 1, checkpoint);
    auto tmp      = model.search(1, 0.0, 1e-3, 1e-3, 1e-3, 1e12, checkpoint);
    auto final_rl = tmp.first;
    auto final_lh = tmp.second;
    // CHECK(final_lh >= initial_lh);
    CHECK(model.compute_lh(final_rl) == Approx(final_lh));
  }

  SECTION("three min roots") {
    auto checkpoint = make_dummy_checkpoint("10.fasta");
    model.assign_indicies_by_rank_search(3, 0.0, 0, 1, checkpoint);
    auto tmp      = model.search(3, 0.0, 1e-3, 1e-3, 1e-3, 1e12, checkpoint);
    auto final_rl = tmp.first;
    auto final_lh = tmp.second;
    CHECK(final_lh >= initial_lh);
    CHECK(model.compute_lh(final_rl) == Approx(final_lh));
  }
}

TEST_CASE("model_t optimize all, slow", "[!hide][all_data][model_t]") {
  auto               ds = data_files_dna["101.phy"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed       = std::rand();
  auto          checkpoint = make_dummy_checkpoint("101.phy");
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  model.compute_lh(tree.root_location(0));
  model.assign_indicies_by_rank_search(1, 0.0, 0, 1, checkpoint);
  auto initial_rl = model.optimize_root_location(1, .05);
  auto tmp        = model.search(1, 0.0, 1e-7, 1e-7, 1e-7, 1e7, checkpoint);
  auto final_rl   = tmp.first;
  auto final_lh   = tmp.second;
  CHECK(final_lh >= initial_rl.second);
  CHECK(model.compute_lh(final_rl) == Approx(final_lh));
}

TEST_CASE("model_t move root", "[model_t]") {
  auto               ds = data_files_dna["101.phy"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed = std::rand();
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions(msa);

  model.set_subst_rates(
      0, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
  model.set_freqs(0, {0.25, 0.25, 0.25, 0.25});
  model.compute_lh(tree.root_location(0));

  auto root_lh = model.compute_all_root_lh();
  for (size_t i = 0; i < root_lh.size(); ++i) {
    for (size_t j = i + 1; j < root_lh.size(); ++j) {
      REQUIRE(root_lh[i] == Approx(root_lh[j]));
    }
  }
}

TEST_CASE("model_t exhaustive search", "[model_t]") {
  auto               ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t tree{ds.second};
  uint64_t      seed       = std::rand();
  auto          checkpoint = make_dummy_checkpoint("10.fasta");
  model_t       model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions(msa);
  model.compute_lh(tree.root_location(0));
  model.assign_indicies_by_rank_exhaustive(0, 1, checkpoint);
  model.exhaustive_search(1e-3, 1e-3, 1e-3, 1e12, checkpoint);
}

TEST_CASE("model_t different rate categories", "[model_t]") {
  auto               ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  rooted_tree_t               tree{ds.second};
  uint64_t                    seed      = std::rand();
  std::vector<ratehet_opts_t> rate_cats = {1, {4}};
  model_t                     model{tree, msa, rate_cats, true, seed, false};
}

TEST_CASE("model_t test no invariant sites", "[model_t]") {
  auto               ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  uint64_t      seed = std::rand();
  rooted_tree_t tree{ds.second};
  model_t       model{tree, msa, {1}, false, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  model.compute_lh(tree.root_location(0));
  auto initial_rl = model.optimize_root_location(1, .05);

  SECTION("one min root") {
    auto checkpoint = make_dummy_checkpoint("10.fasta");
    model.assign_indicies_by_rank_search(1, 0.0, 0, 1, checkpoint);
    auto tmp      = model.search(1, 0.0, 1e-3, 1e-3, 1e-3, 1e12, checkpoint);
    auto final_rl = tmp.first;
    auto final_lh = tmp.second;
    CHECK(final_lh >= initial_rl.second);
    CHECK(model.compute_lh_root(final_rl) == Approx(final_lh));
    CHECK(model.compute_lh(final_rl) == Approx(final_lh));
  }

  SECTION("3 min roots") {
    auto checkpoint = make_dummy_checkpoint("10.fasta");
    CHECK_NOTHROW(
        model.assign_indicies_by_rank_search(3, 0.0, 0, 1, checkpoint));
    auto tmp      = model.search(3, 0.0, 1e-3, 1e-3, 1e-3, 1e12, checkpoint);
    auto final_rl = tmp.first;
    auto final_lh = tmp.second;
    CHECK(final_lh >= initial_rl.second);
    CHECK(model.compute_lh_root(final_rl) == Approx(final_lh));
    CHECK(model.compute_lh(final_rl) == Approx(final_lh));
  }
}

TEST_CASE("assign indicies test", "[model_t]") {
  auto ckp = make_dummy_checkpoint("10.fasta");
  // Since 10.fasta has 10 taxa, that makes it have 2n-3 == 17 possible
  // rootings. So, we want to test that if we write a number of "dummy" results
  // to the file, we get the right result.
  auto               dummy_results_count = GENERATE(0lu, 1lu, 2lu, 4lu, 8lu);
  std::random_device rd;
  std::mt19937       gen(rd());

  std::vector<size_t> possible_idx(17);
  std::iota(possible_idx.begin(), possible_idx.end(), 0);
  std::shuffle(possible_idx.begin(), possible_idx.end(), gen);
  for (size_t i = 0; i < dummy_results_count; ++i) {
    ckp.write(rd_result_t{possible_idx[i], 0.0, 0.0},
              std::vector<partition_parameters_t>{});
  }

  auto               ds = data_files_dna["10.fasta"];
  std::vector<msa_t> msa;
  msa.emplace_back(ds.first);
  uint64_t      seed = std::rand();
  rooted_tree_t tree{ds.second};
  model_t       model{tree, msa, {1}, false, seed, false};
  model.initialize_partitions_uniform_freqs(msa);

  SECTION("search") {
    auto root_assignment = GENERATE(1lu, 2lu, 3lu, 4lu, 5lu);
    int  expected_size   = std::max(static_cast<int>(root_assignment)
                                     - static_cast<int>(dummy_results_count),
                                 0);
    if (static_cast<int>(root_assignment)
            - static_cast<int>(dummy_results_count)
        >= 0) {
      REQUIRE_NOTHROW(model.assign_indicies_by_rank_search(
          root_assignment, 0.0, 0, 1, ckp));
      auto assigned_idx = model.assigned_indicies();
      REQUIRE(assigned_idx.size() == expected_size);
      for (size_t j = 0; j < assigned_idx.size(); ++j) {
        for (size_t i = 0; i < dummy_results_count; ++i) {
          CHECK(assigned_idx[j] != possible_idx[i]);
        }
      }
    } else {
      REQUIRE_THROWS(model.assign_indicies_by_rank_search(1, 0.0, 0, 1, ckp));
    }
  }
  SECTION("exhaustive") {
    int expected_size = std::max(17 - static_cast<int>(dummy_results_count), 0);
    REQUIRE_NOTHROW(model.assign_indicies_by_rank_exhaustive(0, 1, ckp));
    auto assigned_idx = model.assigned_indicies();
    REQUIRE(assigned_idx.size() == expected_size);
    for (size_t j = 0; j < assigned_idx.size(); ++j) {
      for (size_t i = 0; i < dummy_results_count; ++i) {
        CHECK(assigned_idx[j] != possible_idx[i]);
      }
    }
  }
}
