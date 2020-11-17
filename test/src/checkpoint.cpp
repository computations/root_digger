#include "data.hpp"
#include "test_util.hpp"
#include <algorithm>
#include <catch2/catch.hpp>
#include <checkpoint.hpp>
#include <debug.h>
#include <random>
#include <unistd.h>
#include <vector>

std::string make_checkpoint_filename() {
  std::random_device rd;
  uint64_t nonce =
      (static_cast<uint64_t>(rd()) << 32) | static_cast<uint64_t>(rd());
  return std::string("/tmp/checkpoint_test_") + base_58_encode(nonce);
}

checkpoint_t make_and_init_checkpoint() {
  std::string checkpoint_filename = make_checkpoint_filename();
  checkpoint_t ckp(checkpoint_filename);
  cli_options_t cli_options;
  ckp.save_options(cli_options);
  return ckp;
}

TEST_CASE("checkpoint_t constructor", "[checkpoint_t]") {
  std::string checkpoint_filename = make_checkpoint_filename();
  checkpoint_t ckp(checkpoint_filename);
  REQUIRE(access(ckp.get_filename().c_str(), F_OK) != -1);
}

TEST_CASE("checkpoint_t multiple checkpoints", "[checkpoint_t]") {
  std::string checkpoint_filename = make_checkpoint_filename();
  checkpoint_t ckp1(checkpoint_filename);
  REQUIRE(access(ckp1.get_filename().c_str(), F_OK) != -1);
  checkpoint_t ckp2(checkpoint_filename);
  CHECK(ckp2.existing_checkpoint());
}

TEST_CASE("checkpoint_t writing and reading cli_options") {
  std::string checkpoint_filename = make_checkpoint_filename();
  checkpoint_t ckp1(checkpoint_filename);
  REQUIRE(access(ckp1.get_filename().c_str(), F_OK) != -1);
  SECTION("default options") {
    cli_options_t cli_options;
    ckp1.save_options(cli_options);

    checkpoint_t ckp2(checkpoint_filename);
    cli_options_t written_options;
    ckp2.load_options(written_options);
    CHECK(written_options == cli_options);
  }
  SECTION("non-default options") {
    cli_options_t cli_options;
    cli_options.msa_filename = "red roses really like to smell good";
    cli_options.rate_cats = {{1}, {1}, {3}};
    ckp1.save_options(cli_options);

    checkpoint_t ckp2(checkpoint_filename);
    cli_options_t written_options;
    ckp2.load_options(written_options);
    CHECK(written_options == cli_options);
  }
  SECTION("changed options options") {
    cli_options_t cli_options;
    cli_options.msa_filename = "red roses really like to smell good";
    cli_options.rate_cats = {1, 1, 3};
    ckp1.save_options(cli_options);

    cli_options.msa_filename = "this is not the original string";

    checkpoint_t ckp2(checkpoint_filename);
    cli_options_t written_options;
    ckp2.load_options(written_options);
    CHECK(written_options != cli_options);
  }
}

TEST_CASE("checkpoint_t writing and reading results", "[checkpoint_t]") {
  checkpoint_t ckp = make_and_init_checkpoint();
  SECTION("one result") {
    ckp.write(rd_result_t{}, std::vector<partition_parameters_t>{});
    auto results = ckp.read_results();
    CHECK(results.size() == 1);
  }
  SECTION("many results") {
    for (size_t i = 0; i < 1000; ++i) {
      ckp.write(rd_result_t{}, std::vector<partition_parameters_t>{});
    }
    auto results = ckp.read_results();
    CHECK(results.size() == 1000);
  }
}

TEST_CASE("checkpoint_t checking indicies", "[checkpoint_t]") {
  checkpoint_t ckp = make_and_init_checkpoint();
  SECTION("one index") {
    ckp.write(rd_result_t{0, 0.0, 0.0}, std::vector<partition_parameters_t>{});
    auto idx = ckp.completed_indicies();
    REQUIRE(idx.size() == 1);
    CHECK(idx[0] == 0);
  }
  SECTION("generator section") {
    auto total_idx = GENERATE(1lu, 2lu, 4lu, 5lu, 6lu, 7lu, 8lu, 9lu, 10lu);
    for (size_t i = 0; i < total_idx; ++i) {
      ckp.write(rd_result_t{i, 0.0, 0.0},
                std::vector<partition_parameters_t>{});
    }
    auto load_idx = ckp.completed_indicies();
    CHECK(load_idx.size() == total_idx);
    for (size_t i = 0; i < total_idx; ++i) {
      CHECK(std::find(load_idx.begin(), load_idx.end(), i) != load_idx.end());
    }
  }
}
