extern "C" {
#include <libpll/pll.h>
}
#include "data.hpp"
#include <catch2/catch.hpp>
#include <msa.hpp>

TEST_CASE("parse msa", "[msa]") {
  for (auto &&ds : data_files_dna) {
    msa_t msa{ds};
    REQUIRE(msa.states() == 4);
    REQUIRE(msa.map() == pll_map_nt);
  }
}
