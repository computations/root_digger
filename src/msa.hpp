#ifndef __RD_MSA_HPP_
#define __RD_MSA_HPP_

extern "C" {
#include <libpll/pll.h>
}

#include <string>

pll_msa_t *parse_compress_msa_file(const std::string &msa_filename);
pll_msa_t *parse_msa_file(const std::string &msa_filename);

#endif
