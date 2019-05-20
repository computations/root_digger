extern "C" {
#include <libpll/pll.h>
}
#include "data.hpp"
#include <catch2/catch.hpp>
#include <msa.hpp>

TEST_CASE("parse msa", "[msa]") {
  for (auto &&ds : data_files) {
    auto msa = parse_msa_file(ds);
    REQUIRE(msa != nullptr);
    pll_msa_destroy(msa);
  }
}
