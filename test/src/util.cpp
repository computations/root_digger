#include "data.hpp"
#include "test_util.hpp"
#include <catch2/catch.hpp>
#include <random>
#include <util.hpp>

TEST_CASE("cli_options_t comparison operators", "[cli_options_t]"){
  cli_options_t cli1;
  cli1.msa_filename = "what in the world is a kangaroo doing in the room";
  cli1.seed = 12312;
  cli_options_t cli2;
  cli2.msa_filename = "what in the world is a kangaroo doing in the room";
  cli2.seed = 12312;
  SECTION("operator=="){
    CHECK(cli1 == cli2);
  }
  SECTION("operator!="){
    cli1.msa_filename = "this is now a different string than it was origionally";
    CHECK(cli1 != cli2);
  }
}
