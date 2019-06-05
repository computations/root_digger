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
    for (size_t i = 0; i < tree.root_count(); ++i) {
      CHECK(tree.root_location(i).edge->length > 0.0);
    }
  }
}

TEST_CASE("rooted_tree string constructor no file", "[rooted_tree_t]") {
  REQUIRE_THROWS(rooted_tree_t{"not_a_tree_file"});
}

TEST_CASE("rooted_tree copy constructor", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2{tree1};
    CHECK(tree1.branch_count() == tree2.branch_count());
    CHECK(tree1.inner_count() == tree2.inner_count());
    CHECK(tree1.tip_count() == tree2.tip_count());
  }
}

TEST_CASE("rooted_tree copy asignment constructor", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2 = tree1;
    CHECK(tree1.branch_count() == tree2.branch_count());
    CHECK(tree1.inner_count() == tree2.inner_count());
    CHECK(tree1.tip_count() == tree2.tip_count());
  }
}

TEST_CASE("rooted_tree move asignment constructor", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2 = std::move(tree1);
    CHECK(tree2.branch_count() > 0);
    CHECK(tree2.inner_count() > 0);
    CHECK(tree2.tip_count() > 0);
  }
}

TEST_CASE("rooted_tree_t label map", "[rooted_tree_t]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree{ds.second};
    auto lm = tree.label_map();
    CHECK(lm.size() == tree.tip_count());
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
      CHECK(ops.size() > 0);
      CHECK(pmatrices.size() > 0);
      CHECK(branches.size() > 0);
      for (auto brlen : branches) {
        CHECK(brlen > 0.0);
      }
    }
  }
}

TEST_CASE("rooted_tree_t generate operations, known tree",
          "[rooted_tree_t][regression]") {
  rooted_tree_t tree{data_files_dna[0].second};
  std::vector<pll_operation_t> ops;
  std::vector<unsigned int> pmatrices;
  std::vector<double> branches;

  CHECK(std::string(tree.root_location(3).edge->label) == "n2");

  GENERATE_AND_UNPACK_OPS(tree, tree.root_location(3), ops, pmatrices,
                          branches);
  CHECK(ops.size() == 3);

  CHECK(ops[0].parent_clv_index == 4);
  CHECK(ops[0].parent_scaler_index == 0);
  CHECK(ops[0].child1_clv_index == 0);
  CHECK(ops[0].child1_scaler_index == -1);
  CHECK(ops[0].child1_matrix_index == 0);
  CHECK(ops[0].child2_clv_index == 1);
  CHECK(ops[0].child2_scaler_index == -1);
  CHECK(ops[0].child2_matrix_index == 1);

  CHECK(ops[1].parent_clv_index == 5);
  CHECK(ops[1].parent_scaler_index == 1);
  CHECK(ops[1].child1_clv_index == 2);
  CHECK(ops[1].child1_scaler_index == -1);
  CHECK(ops[1].child1_matrix_index == 2);
  CHECK(ops[1].child2_clv_index == 3);
  CHECK(ops[1].child2_scaler_index == -1);
  CHECK(ops[1].child2_matrix_index == 3);

  CHECK(ops[2].parent_clv_index == 6);
  CHECK(ops[2].parent_scaler_index == 2);
  CHECK(ops[2].child1_clv_index == 4);
  CHECK(ops[2].child1_scaler_index == 0);
  CHECK(ops[2].child1_matrix_index == 4);
  CHECK(ops[2].child2_clv_index == 5);
  CHECK(ops[2].child2_scaler_index == 1);
  CHECK(ops[2].child2_matrix_index == 5);
}

TEST_CASE("rooted_tree_t root operations", "[rooted_tree_t][root_by]") {
  for (auto &ds : data_files_dna) {
    rooted_tree_t tree{ds.second};
    for (size_t i = 0; i < tree.root_count(); ++i) {
      tree.root_by(tree.root_location(i));
      tree.unroot();
    }
  }
}

TEST_CASE("rooted_tree_t newick", "[rooted_tree_t]") {
  rooted_tree_t tree{data_files_dna[0].second};
  REQUIRE(tree.root_count() == 5);
  auto rl1 = tree.root_location(0);
  rl1.brlen_ratio = 0.25;
  tree.root_by(rl1);
  CHECK("(b:0.025000,((c:0.100000,d:0.100000)n2:0.550000,a:0.100000)n1:0."
        "075000):0.0;" == tree.newick());

  rl1.brlen_ratio = 0.75;
  tree.root_by(rl1);
  CHECK("(b:0.075000,((c:0.100000,d:0.100000)n2:0.550000,a:0.100000)n1:0."
        "025000):0.0;" == tree.newick());

  auto rl2 = tree.root_location(1);
  rl2.brlen_ratio = 0.25;
  tree.root_by(rl2);
  CHECK("(a:0.025000,(b:0.100000,(c:0.100000,d:0.100000)n2:0.550000)n1:0."
        "075000):0.0;" == tree.newick());

  rl2.brlen_ratio = 0.75;
  tree.root_by(rl2);
  CHECK("(a:0.075000,(b:0.100000,(c:0.100000,d:0.100000)n2:0.550000)n1:0."
        "025000):0.0;" == tree.newick());

  auto rl3 = tree.root_location(2);
  rl3.brlen_ratio = 0.25;
  tree.root_by(rl3);
  CHECK("((c:0.100000,d:0.100000)n2:0.137500,(a:0.100000,b:0.100000)n1:0."
        "412500):0.0;" == tree.newick());

  rl3.brlen_ratio = 0.75;
  tree.root_by(rl3);
  CHECK("((c:0.100000,d:0.100000)n2:0.412500,(a:0.100000,b:0.100000)n1:0."
        "137500):0.0;" == tree.newick());

  auto rl4 = tree.root_location(3);
  rl4.brlen_ratio = 0.25;
  tree.root_by(rl4);
  CHECK("(c:0.025000,(d:0.100000,(a:0.100000,b:0.100000)n1:0.550000)n2:0."
        "075000):0.0;" == tree.newick());

  rl4.brlen_ratio = 0.75;
  tree.root_by(rl4);
  CHECK("(c:0.075000,(d:0.100000,(a:0.100000,b:0.100000)n1:0.550000)n2:0."
        "025000):0.0;" == tree.newick());

  auto rl5 = tree.root_location(4);
  rl5.brlen_ratio = 0.25;
  tree.root_by(rl5);
  CHECK("(d:0.025000,((a:0.100000,b:0.100000)n1:0.550000,c:0.100000)n2:0."
        "075000):0.0;" == tree.newick());

  rl5.brlen_ratio = 0.75;
  tree.root_by(rl5);
  CHECK("(d:0.075000,((a:0.100000,b:0.100000)n1:0.550000,c:0.100000)n2:0."
        "025000):0.0;" == tree.newick());
}
