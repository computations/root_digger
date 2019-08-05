#include "data.hpp"

#define __stringify(s) std::string(#s)
#define __stringify_base(s) __stringify(s)

std::vector<std::string> data_keys = {
    "single",
    "101.phy",
};

std::map<std::string, std::pair<std::string, std::string>>
    data_files_dna = {
        {"single",
         {__stringify_base(DATA_DIRECTORY_DNA_ABS) + "single.phy",
          __stringify_base(DATA_DIRECTORY_TREE_ABS) + "single.tree"}},
        {"101.phy",
         {__stringify_base(DATA_DIRECTORY_DNA_ABS) + "101.phy",
          __stringify_base(DATA_DIRECTORY_TREE_ABS) + "101.tree"}},
};

std::map<std::string, std::string> check_trees = {
    {"sanity_check1",
     __stringify_base(DATA_DIRECTORY_TREE_ABS) + "sanity_check1.tree"},
    {"sanity_check2",
     __stringify_base(DATA_DIRECTORY_TREE_ABS) + "sanity_check2.tree"},
    {"sanity_check3",
     __stringify_base(DATA_DIRECTORY_TREE_ABS) + "sanity_check3.tree"},
};
