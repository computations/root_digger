extern "C" {
#include <libpll/pll.h>
}
#include "data.hpp"
#include <catch2/catch.hpp>
#include <msa.hpp>

TEST_CASE("msa_t parse msa", "[msa_t]") {
  for (auto &kv : data_files_dna) {
    auto &ds = kv.second;
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
  SECTION("Model strings") {
    SECTION("Frequency only") {
      SECTION("emperical") {
        std::string line{"DNA+F,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+F" == pi.model_name);
        CHECK("DNA" == pi.model.subst_str);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.freq_opts.type == param_type::emperical);
      }
      SECTION("estimate") {
        std::string line{"DNA+FO,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+FO" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.freq_opts.type == param_type::estimate);
      }
      SECTION("equal") {
        std::string line{"DNA+FE,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+FE" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.freq_opts.type == param_type::equal);
      }
      SECTION("user") {
        SECTION("cli") {
          std::string line{
              "DNA+FU{0.25/0.25/0.25/0.25},PART_0=123-4123, 5122-12411"};
          auto pi = parse_partition_info(line);
          CHECK("DNA+FU{0.25/0.25/0.25/0.25}" == pi.model_name);
          CHECK("PART_0" == pi.partition_name);
          CHECK(123 == pi.parts[0].first);
          CHECK(4123 == pi.parts[0].second);
          CHECK(5122 == pi.parts[1].first);
          CHECK(12411 == pi.parts[1].second);
          CHECK(pi.model.freq_opts.type == param_type::user);
        }
      }
    }
    SECTION("Invar Only"){
      SECTION("estimate"){
        std::string line{"DNA+I,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+I" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.invar_opts.type == param_type::estimate);
      }
      SECTION("emperical"){
        std::string line{"DNA+IC,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+IC" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.invar_opts.type == param_type::emperical);
      }
      SECTION("user"){
        std::string line{"DNA+IU{0.25},PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+IU{0.25}" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.invar_opts.type == param_type::user);
        CHECK(pi.model.invar_opts.user_prop == 0.25);
      }
    }
    SECTION("RateHet Only"){
      SECTION("emperical"){
        std::string line{"DNA+G,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+G" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.ratehet_opts.type == param_type::estimate);
        CHECK(pi.model.ratehet_opts.rate_cats == 4);
      }
      SECTION("emperical, with n"){
        std::string line{"DNA+G2,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+G2" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.ratehet_opts.type == param_type::estimate);
        CHECK(pi.model.ratehet_opts.rate_cats == 2);
      }
      SECTION("user"){
        std::string line{"DNA+G2{0.25},PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+G2{0.25}" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.ratehet_opts.type == param_type::user);
        CHECK(pi.model.ratehet_opts.rate_cats == 2);
        CHECK(pi.model.ratehet_opts.alpha == 0.25);
      }
      SECTION("median"){
        std::string line{"DNA+GA,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+GA" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.ratehet_opts.type == param_type::estimate);
        CHECK(pi.model.ratehet_opts.rate_cats == 4);
      }
    }
    SECTION("All together"){
        std::string line{"DNA+G2{0.25}+F+I,PART_0=123-4123, 5122-12411"};
        auto pi = parse_partition_info(line);
        CHECK("DNA+G2{0.25}+F+I" == pi.model_name);
        CHECK("PART_0" == pi.partition_name);
        CHECK(123 == pi.parts[0].first);
        CHECK(4123 == pi.parts[0].second);
        CHECK(5122 == pi.parts[1].first);
        CHECK(12411 == pi.parts[1].second);
        CHECK(pi.model.ratehet_opts.type == param_type::user);
        CHECK(pi.model.ratehet_opts.rate_cats == 2);
        CHECK(pi.model.ratehet_opts.alpha == 0.25);
        CHECK(pi.model.invar_opts.type == param_type::estimate);
        CHECK(pi.model.freq_opts.type == param_type::emperical);
    }
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

TEST_CASE("msa_t partition datafile", "[msa_t]") {
  auto ds = data_files_dna["101.phy"];
  SECTION("single partition, one range") {
    msa_t msa{ds.first};
    msa_partitions_t parts{parse_partition_info("DNA, PART_0 = 0-100")};
    auto parted_msa = msa.partition(parts);
    CHECK(parted_msa.size() == parts.size());
    CHECK(parted_msa[0].length() == 101);
  }
  SECTION("single partition, two ranges") {
    msa_t msa{ds.first};
    msa_partitions_t parts{
        parse_partition_info("DNA, PART_0 = 0-100, 200-300")};
    auto parted_msa = msa.partition(parts);
    CHECK(parted_msa.size() == parts.size());
    CHECK(parted_msa[0].length() == 202);
  }
  SECTION("multiple partitions, one range") {
    msa_t msa{ds.first};
    msa_partitions_t parts{parse_partition_info("DNA, PART_0 = 0-100"),
                           parse_partition_info("DNA, PART_1 = 200-300")};
    auto parted_msa = msa.partition(parts);
    CHECK(parted_msa.size() == parts.size());
    CHECK(parted_msa[0].length() == 101);
    CHECK(parted_msa[1].length() == 101);
  }
  SECTION("multiple partitions, multiple ranges") {
    msa_t msa{ds.first};
    msa_partitions_t parts{
        parse_partition_info("DNA, PART_0 = 0-100, 500-520"),
        parse_partition_info("DNA, PART_1 = 200-300, 400-500")};
    auto parted_msa = msa.partition(parts);
    CHECK(parted_msa.size() == parts.size());
    CHECK(parted_msa[0].length() == 122);
    CHECK(parted_msa[1].length() == 202);
  }
}
