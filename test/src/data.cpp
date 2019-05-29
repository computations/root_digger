#include "data.hpp"

#define __stringify(s) std::string(#s)
#define __stringify_base(s) __stringify(s)

std::vector<std::pair<std::string, std::string>> data_files_dna = {
    {__stringify_base(DATA_DIRECTORY_DNA_ABS) + "single.phy",
     __stringify_base(DATA_DIRECTORY_TREE_ABS) + "single.tree"},
    {__stringify_base(DATA_DIRECTORY_DNA_ABS) + "101.phy",
     __stringify_base(DATA_DIRECTORY_TREE_ABS) + "101.tree"},
};
