#ifndef __RD_TEST_DATA_HPP_
#define __RD_TEST_DATA_HPP_
#include <vector>
#include <string>

#define __stringify(s) std::string(#s)
#define __stringify_base(s) __stringify(s)

std::vector<std::string> data_files = {
    __stringify_base(DATA_DIRECTORY_ABS) + "101.phy",
};

#endif
