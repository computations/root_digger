extern "C" {
#include <libpll/pll.h>
}
#include "data.hpp"
#include <catch2/catch.hpp>
#include <msa.hpp>

TEST_CASE("msa_t parse msa", "[msa_t]") {
  for (auto &&ds : data_files_dna) {
    msa_t msa{ds.first};
    REQUIRE(msa.states() == 4);
    REQUIRE(msa.map() == pll_map_nt);
  }
}

TEST_CASE("msa_t parse partition line", "[msa_t]") {
  SECTION("no errors") {
    std::string line{"DNA, PART_0 = 123-4123"};
    auto pi = parse_partition_info(line);
    CHECK("DNA" == pi.model_name);
    CHECK("PART_0" == pi.partition_name);
    CHECK(123 == pi.parts[0].first);
    CHECK(4123 == pi.parts[0].second);
  }
  SECTION("no spaces") {
    std::string line{"DNA,PART_0=123-4123"};
    auto pi = parse_partition_info(line);
    CHECK("DNA" == pi.model_name);
    CHECK("PART_0" == pi.partition_name);
    CHECK(123 == pi.parts[0].first);
    CHECK(4123 == pi.parts[0].second);
  }
  SECTION("multiple ranges") {
    std::string line{"DNA,PART_0=123-4123, 5122-12411"};
    auto pi = parse_partition_info(line);
    CHECK("DNA" == pi.model_name);
    CHECK("PART_0" == pi.partition_name);
    CHECK(123 == pi.parts[0].first);
    CHECK(4123 == pi.parts[0].second);
    CHECK(5122 == pi.parts[1].first);
    CHECK(12411 == pi.parts[1].second);
  }
  SECTION("error: missing comma") {
    std::string line{"DNA PART_0 = 123-4123"};
    CHECK_THROWS(parse_partition_info(line));
  }
  SECTION("error: missing =") {
    std::string line{"DNA, PART_0  123-4123"};
    CHECK_THROWS(parse_partition_info(line));
  }
  SECTION("error: missing -") {
    std::string line{"DNA, PART_0 = 1234123"};
    CHECK_THROWS(parse_partition_info(line));
  }
  SECTION("error: missing -") {
    std::string line{"DNA, PART_0 = 123=4123"};
    CHECK_THROWS(parse_partition_info(line));
  }
  SECTION("error: missing model name") {
    std::string line{", PART_0 = 123-4123"};
    CHECK_THROWS(parse_partition_info(line));
  }
}
