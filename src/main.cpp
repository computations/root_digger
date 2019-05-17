extern "C" {
#include <libpll/pll.h>
}
#include <getopt.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <vector>

typedef std::vector<double> model_params_t;

std::string read_file_contents(std::ifstream& infile){
  std::string str;
  infile.seekg(0, std::ios::end);
  str.resize(infile.tellg());
  infile.seekg(0, std::ios::beg);
  infile.read(&str[0], str.size());
  return str;
}

double parse_param(std::string::const_iterator begin,
                   std::string::const_iterator end) {
  try {
    std::string tmp{begin, end};
    return std::stod(tmp);
  } catch (const std::invalid_argument &) {
    throw std::invalid_argument("Failed to parse model parameter string");
  } catch (const std::out_of_range &) {
    throw std::out_of_range(
        "Model parameter string contains a value which is too large");
  }
}

model_params_t parse_model_params(const std::string &model_string) {
  model_params_t model_params;
  auto substr_begin = model_string.begin();

  for (auto substr_end = model_string.begin();
       substr_begin != model_string.end(); ++substr_end) {

    if (*substr_end == ',') {
      model_params.push_back(parse_param(substr_begin, substr_end));
      substr_begin = substr_end + 1;
    } else if (*substr_end == '\0') {
      model_params.push_back(parse_param(substr_begin, substr_end));
      break;
    }
  }

  return model_params;
}

model_params_t parse_model_file(const std::string& model_filename){
  std::ifstream model_file(model_filename);

  return parse_model_params(read_file_contents(model_file));
}

int main(int argv, char **argc) {
  static struct option long_opts[] = {
      {"msa", required_argument, 0, 0},
      {"tree", required_argument, 0, 0},
      {"model", required_argument, 0, 0},
      {0, 0, 0, 0},
  };

  char c;
  int index = 0;
  std::string msa_filename;
  std::string tree_filename;
  std::string model_filename;
  while ((c = getopt_long_only(argv, argc, "", long_opts, &index))) {
    switch (index){
      case 0:
        msa_filename = optarg;
        break;
      case 1:
        tree_filename = optarg;
        break;
      case 2:
        model_filename = optarg;
        break;
    }
  }

  return 0;
}
