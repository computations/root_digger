extern "C" {
#include <libpll/pll.h>
}
#include "data.hpp"
#include <catch2/catch.hpp>
#include <debug.h>
#include <tree.hpp>
#include <utility>

root_location_t find_node(const rooted_tree_t &tree, const std::string &label) {
  for (size_t i = 0; i < tree.root_count(); ++i) {
    if (tree.root_location(i).edge->label == label)
      return tree.root_location(i);
  }
  throw std::invalid_argument("Could not find the correct node");
}

TEST_CASE("rooted_tree_t string constructor", "[rooted_tree_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
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
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2{tree1};
    CHECK(tree1.branch_count() == tree2.branch_count());
    CHECK(tree1.inner_count() == tree2.inner_count());
    CHECK(tree1.tip_count() == tree2.tip_count());
  }
}

TEST_CASE("rooted_tree copy asignment constructor", "[rooted_tree_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2 = tree1;
    CHECK(tree1.branch_count() == tree2.branch_count());
    CHECK(tree1.inner_count() == tree2.inner_count());
    CHECK(tree1.tip_count() == tree2.tip_count());
  }
}

TEST_CASE("rooted_tree move asignment constructor", "[rooted_tree_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    rooted_tree_t tree1{ds.second};
    rooted_tree_t tree2 = std::move(tree1);
    CHECK(tree2.branch_count() > 0);
    CHECK(tree2.inner_count() > 0);
    CHECK(tree2.tip_count() > 0);
  }
}

TEST_CASE("rooted_tree_t label map", "[rooted_tree_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    rooted_tree_t tree{ds.second};
    auto lm = tree.label_map();
    CHECK(lm.size() == tree.tip_count());
    for (const auto &kv : lm) {
      CHECK(kv.second <= tree.tip_count());
    }
  }
}

TEST_CASE("rooted_tree_t generate operations", "[rooted_tree_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
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
  rooted_tree_t tree{data_files_dna["single"].second};
  std::vector<pll_operation_t> ops;
  std::vector<unsigned int> pmatrices;
  std::vector<double> branches;

  auto root_location = find_node(tree, "n2");

  GENERATE_AND_UNPACK_OPS(tree, root_location, ops, pmatrices, branches);
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

TEST_CASE("rooted_tree_t generate root operations", "[rooted_tree_t]") {
  rooted_tree_t tree{data_files_dna["single"].second};
  pll_operation_t op;
  std::vector<unsigned int> pmatrices;
  std::vector<double> branches;

  auto root_location = find_node(tree, "n2");
  {
    auto result = tree.generate_derivative_operations(root_location);
    op = std::move(std::get<0>(result));
    pmatrices = std::move(std::get<1>(result));
    branches = std::move(std::get<2>(result));
  }

  CHECK(op.parent_clv_index == 6);
  CHECK(op.parent_scaler_index == 2);
  CHECK(op.child1_clv_index == 4);
  CHECK(op.child1_scaler_index == 0);
  CHECK(op.child1_matrix_index == 4);
  CHECK(op.child2_clv_index == 5);
  CHECK(op.child2_scaler_index == 1);
  CHECK(op.child2_matrix_index == 5);

  CHECK(pmatrices.size() == 2);
  CHECK(pmatrices[0] == 4);
  CHECK(pmatrices[1] == 5);

  CHECK(branches.size() == 2);
  CHECK(branches[0] == 0.275);
  CHECK(branches[1] == 0.275);
}

TEST_CASE("rooted_tree_t root operations", "[rooted_tree_t][root_by]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
    rooted_tree_t tree{ds.second};
    for (size_t i = 0; i < tree.root_count(); ++i) {
      tree.root_by(tree.root_location(i));
      tree.unroot();
    }
  }
}

TEST_CASE("rooted_tree_t newick", "[rooted_tree_t]") {
  rooted_tree_t tree{data_files_dna["single"].second};
  REQUIRE(tree.root_count() == 5);
  auto rl1 = find_node(tree, "b");
  rl1.brlen_ratio = 0.25;
  tree.root_by(rl1);
  CHECK("(b:0.025000,((c:0.100000,d:0.100000)n2:0.550000,a:0.100000)n1:0."
        "075000):0.0;" == tree.newick());

  rl1.brlen_ratio = 0.75;
  tree.root_by(rl1);
  CHECK("(b:0.075000,((c:0.100000,d:0.100000)n2:0.550000,a:0.100000)n1:0."
        "025000):0.0;" == tree.newick());

  auto rl2 = find_node(tree, "a");
  rl2.brlen_ratio = 0.25;
  tree.root_by(rl2);
  CHECK("(a:0.025000,(b:0.100000,(c:0.100000,d:0.100000)n2:0.550000)n1:0."
        "075000):0.0;" == tree.newick());

  rl2.brlen_ratio = 0.75;
  tree.root_by(rl2);
  CHECK("(a:0.075000,(b:0.100000,(c:0.100000,d:0.100000)n2:0.550000)n1:0."
        "025000):0.0;" == tree.newick());

  auto rl3 = find_node(tree, "n2");
  rl3.brlen_ratio = 0.25;
  tree.root_by(rl3);
  CHECK("((c:0.100000,d:0.100000)n2:0.137500,(a:0.100000,b:0.100000)n1:0."
        "412500):0.0;" == tree.newick());

  rl3.brlen_ratio = 0.75;
  tree.root_by(rl3);
  CHECK("((c:0.100000,d:0.100000)n2:0.412500,(a:0.100000,b:0.100000)n1:0."
        "137500):0.0;" == tree.newick());

  auto rl4 = find_node(tree, "c");
  rl4.brlen_ratio = 0.25;
  tree.root_by(rl4);
  CHECK("(c:0.025000,(d:0.100000,(a:0.100000,b:0.100000)n1:0.550000)n2:0."
        "075000):0.0;" == tree.newick());

  rl4.brlen_ratio = 0.75;
  tree.root_by(rl4);
  CHECK("(c:0.075000,(d:0.100000,(a:0.100000,b:0.100000)n1:0.550000)n2:0."
        "025000):0.0;" == tree.newick());

  auto rl5 = find_node(tree, "d");
  rl5.brlen_ratio = 0.25;
  tree.root_by(rl5);
  CHECK("(d:0.025000,((a:0.100000,b:0.100000)n1:0.550000,c:0.100000)n2:0."
        "075000):0.0;" == tree.newick());

  rl5.brlen_ratio = 0.75;
  tree.root_by(rl5);
  CHECK("(d:0.075000,((a:0.100000,b:0.100000)n1:0.550000,c:0.100000)n2:0."
        "025000):0.0;" == tree.newick());
}

TEST_CASE("rooted_tree_t check move root operations", "[rooted_tree_t]") {
  rooted_tree_t tree{data_files_dna["single"].second};
}

TEST_CASE("rooted_tree_t check derivative vs regular operations",
          "[rooted_tree_t]") {
  rooted_tree_t tree{data_files_dna["single"].second};
  for (size_t i = 0; i < tree.root_count(); ++i) {
    auto root_location = tree.root_location(i);
    std::vector<pll_operation_t> ops_regular;
    std::vector<unsigned int> pmatrices_regular;
    std::vector<double> branches_regular;

    GENERATE_AND_UNPACK_OPS(tree, root_location, ops_regular, pmatrices_regular,
                            branches_regular);

    pll_operation_t op_root;
    std::vector<unsigned int> pmatrices_root;
    std::vector<double> branches_root;

    {
      auto result = tree.generate_derivative_operations(root_location);
      op_root = std::move(std::get<0>(result));
      pmatrices_root = std::move(std::get<1>(result));
      branches_root = std::move(std::get<2>(result));
    }

    auto regular_end_op = *(ops_regular.end() - 1);
    CHECK(op_root.parent_clv_index == regular_end_op.parent_clv_index);
    CHECK(op_root.parent_scaler_index == regular_end_op.parent_scaler_index);

    CHECK(op_root.child1_clv_index == regular_end_op.child1_clv_index);
    CHECK(op_root.child2_clv_index == regular_end_op.child2_clv_index);

    CHECK(op_root.child1_scaler_index == regular_end_op.child1_scaler_index);
    CHECK(op_root.child2_scaler_index == regular_end_op.child2_scaler_index);

    CHECK(op_root.child1_matrix_index == regular_end_op.child1_matrix_index);
    CHECK(op_root.child2_matrix_index == regular_end_op.child2_matrix_index);
  }
}

TEST_CASE("rooted_tree_t sanity check", "[rooted_tree_t]") {
  rooted_tree_t t1(check_trees["sanity_check1"]);
  CHECK(!t1.sanity_check());

  rooted_tree_t t2(check_trees["sanity_check2"]);
  CHECK(!t2.sanity_check());

  rooted_tree_t t3(check_trees["sanity_check3"]);
  CHECK(t3.sanity_check());
}

TEST_CASE("rooted_tree_t annotations, basic", "[rooted_tree_t]") {
  rooted_tree_t t1(data_files_dna["10.fasta"].second);
  for (auto rl : t1.roots()) {
    t1.annotate_branch(rl, "foo", "bar");
    t1.annotate_branch(rl, "fizz", "buzz");
  }
  CHECK(
      t1.newick() ==
      "(((j:0.854700[&&NHX:foo=bar:fizz=buzz],((h:0.983500[&&NHX:foo=bar:fizz="
      "buzz],a:0.224900[&&NHX:foo=bar:fizz=buzz]):0.416200[&&NHX:foo=bar:fizz="
      "buzz],(c:0.540900[&&NHX:foo=bar:fizz=buzz],f:0.422200[&&NHX:foo=bar:"
      "fizz=buzz]):0.785300[&&NHX:foo=bar:fizz=buzz]):0.614100[&&NHX:foo=bar:"
      "fizz=buzz]):0.446100[&&NHX:foo=bar:fizz=buzz],g:0.487400[&&NHX:foo=bar:"
      "fizz=buzz]):0.825200[&&NHX:foo=bar:fizz=buzz],((i:0.569700[&&NHX:foo="
      "bar:fizz=buzz],e:0.366600[&&NHX:foo=bar:fizz=buzz]):0.602800[&&NHX:foo="
      "bar:fizz=buzz],b:0.445900[&&NHX:foo=bar:fizz=buzz]):0.099300[&&NHX:foo="
      "bar:fizz=buzz],d:0.639600[&&NHX:foo=bar:fizz=buzz]):0.0;");
}

TEST_CASE("rooted_tree_t annotations, moving root", "[rooted_tree_t]") {
  rooted_tree_t t1(data_files_dna["10.fasta"].second);
  t1.root_by(t1.roots()[0]);
  for (auto rl : t1.roots()) {
    t1.annotate_branch(rl, "foo", "bar");
    t1.annotate_branch(rl, "fizz", "buzz");
  }
  t1.root_by(t1.roots()[1]);
  t1.unroot();
  CHECK(
      t1.newick() ==
      "((d:0.639600[&&NHX:foo=bar:fizz=buzz],((j:0.854700[&&NHX:foo=bar:fizz="
      "buzz],((h:0.983500[&&NHX:foo=bar:fizz=buzz],a:0.224900[&&NHX:foo=bar:"
      "fizz=buzz]):0.416200[&&NHX:foo=bar:fizz=buzz],(c:0.540900[&&NHX:foo=bar:"
      "fizz=buzz],f:0.422200[&&NHX:foo=bar:fizz=buzz]):0.785300[&&NHX:foo=bar:"
      "fizz=buzz]):0.614100[&&NHX:foo=bar:fizz=buzz]):0.446100[&&NHX:foo=bar:"
      "fizz=buzz],g:0.487400[&&NHX:foo=bar:fizz=buzz]):0.825200[&&NHX:foo=bar:"
      "fizz=buzz]):0.099300[&&NHX:foo=bar:fizz=buzz],(i:0.569700[&&NHX:foo=bar:"
      "fizz=buzz],e:0.366600[&&NHX:foo=bar:fizz=buzz]):0.602800[&&NHX:foo=bar:"
      "fizz=buzz],b:0.445900[&&NHX:foo=bar:fizz=buzz]):0.0;");
}

TEST_CASE("rooted_tree_t annotations, all roots", "[rooted_tree_t]") {
  rooted_tree_t t1(data_files_dna["10.fasta"].second);
  for (auto rl : t1.roots()) {
    t1.root_by(rl);
    t1.annotate_branch(rl, "foo", "bar");
    t1.annotate_branch(rl, "fizz", "buzz");
  }
  t1.unroot();
  CHECK(
      t1.newick() ==
      "((((i:0.569700[&&NHX:foo=bar:fizz=buzz],e:0.366600[&&NHX:foo=bar:fizz="
      "buzz]):0.602800[&&NHX:foo=bar:fizz=buzz],b:0.445900[&&NHX:foo=bar:fizz="
      "buzz]):0.099300[&&NHX:foo=bar:fizz=buzz],d:0.639600[&&NHX:foo=bar:fizz="
      "buzz]):0.825200[&&NHX:foo=bar:fizz=buzz],(j:0.854700[&&NHX:foo=bar:fizz="
      "buzz],((h:0.983500[&&NHX:foo=bar:fizz=buzz],a:0.224900[&&NHX:foo=bar:"
      "fizz=buzz]):0.416200[&&NHX:foo=bar:fizz=buzz],(c:0.540900[&&NHX:foo=bar:"
      "fizz=buzz],f:0.422200[&&NHX:foo=bar:fizz=buzz]):0.785300[&&NHX:foo=bar:"
      "fizz=buzz]):0.614100[&&NHX:foo=bar:fizz=buzz]):0.446100[&&NHX:foo=bar:"
      "fizz=buzz],g:0.487400[&&NHX:foo=bar:fizz=buzz]):0.0;");
}
