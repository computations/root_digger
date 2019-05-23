extern "C" {
#include <libpll/pll.h>
}
#include "data.hpp"
#include <catch2/catch.hpp>
#include <tree.hpp>
#include <utility>

TEST_CASE("rooted_tree_t string constructor", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree{ds.second};
    CHECK(tree.root_count() > 0);
  }
}

TEST_CASE("rooted_tree string constructor no file", "[rooted_tree_t]") {
  REQUIRE_THROWS(rooted_tree_t{"not_a_tree_file"});
}

TEST_CASE("rooted_tree copy constructor", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2{tree1};
    CHECK(tree1.branches() == tree2.branches());
    CHECK(tree1.inner_count() == tree2.inner_count());
    CHECK(tree1.tip_count() == tree2.tip_count());
  }
}

TEST_CASE("rooted_tree copy asignment constructor", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2 = tree1;
    CHECK(tree1.branches() == tree2.branches());
    CHECK(tree1.inner_count() == tree2.inner_count());
    CHECK(tree1.tip_count() == tree2.tip_count());
  }
}

TEST_CASE("rooted_tree move asignment constructor", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2 = std::move(tree1);
    CHECK(tree2.branches() > 0);
    CHECK(tree2.inner_count() > 0);
    CHECK(tree2.tip_count() > 0);
  }
}

TEST_CASE("rooted_tree_t label map", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree{ds.second};
    auto lm = tree.label_map();
    for (const auto &kv : lm) {
      CHECK(kv.second <= tree.tip_count());
    }
  }
}

TEST_CASE("rooted_tree_t generate operations", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree{ds.second};
    std::vector<pll_operation_t> ops;
    std::vector<unsigned int> pmatrices;
    std::vector<double> branches;

    for (size_t i = 0; i < tree.root_count(); ++i) {
      GENERATE_AND_UNPACK_OPS(tree, tree.root_location(i), ops, pmatrices,
                              branches);
    }
  }
}
