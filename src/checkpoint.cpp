#include "checkpoint.hpp"
#include "debug.h"
#include "util.hpp"
#include <fcntl.h>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <functional>

template <> size_t write(int fd, const std::string &str) {
  size_t string_size  = str.size();
  size_t total_writen = 0;
  auto   result       = write(fd, &string_size, sizeof(string_size));
  if (result < static_cast<ssize_t>(sizeof(string_size))) {
    throw checkpoint_write_failure{
        "Failed to write a string to the checkpoint file"};
  }
  total_writen += static_cast<size_t>(result);

  result = write(fd, str.data(), string_size * sizeof(std::string::value_type));
  total_writen += static_cast<size_t>(result);

  return total_writen;
}

template <> size_t read(int fd, std::string &str) {
  size_t total_read;
  size_t string_size = 0;

  total_read = read(fd, string_size);

  if (string_size == 0) { return total_read; }

  char *tmp_str = (char *)calloc(string_size + 1, sizeof(char));
  auto  result  = read(fd, tmp_str, string_size * sizeof(char));
  if (result < 0) { throw checkpoint_read_failure{"Failed to read a string"}; }
  total_read += static_cast<size_t>(result);

  std::string tmp(tmp_str);
  std::swap(tmp, str);
  free(tmp_str);

  return total_read;
}

template <> size_t write(int fd, const cli_options_t &options) {
  size_t total_written = 0;
  total_written += write(fd, options.msa_filename);
  total_written += write(fd, options.tree_filename);
  total_written += write(fd, options.prefix);
  total_written += write(fd, options.model_filename);
  total_written += write(fd, options.freqs_filename);
  total_written += write(fd, options.partition_filename);
  total_written += write(fd, options.data_type);
  total_written += write(fd, options.model_string);
  total_written += write(fd, options.rate_cats);
  total_written += write(fd, options.rate_category_types);
  total_written += write(fd, options.seed);
  total_written += write(fd, options.min_roots);
  total_written += write(fd, options.threads);
  total_written += write(fd, options.root_ratio);
  total_written += write(fd, options.abs_tolerance);
  total_written += write(fd, options.factor);
  total_written += write(fd, options.br_tolerance);
  total_written += write(fd, options.bfgs_tol);
  // total_written += write(fd, options.states);
  total_written += write(fd, options.silent);
  total_written += write(fd, options.exhaustive);
  total_written += write(fd, options.echo);
  total_written += write(fd, options.invariant_sites);
  total_written += write(fd, options.early_stop);
  debug_print(
      EMIT_LEVEL_DEBUG, "Wrote %lu bytes for cli_options", total_written);
  return total_written;
}

template <> size_t read(int fd, cli_options_t &options) {
  size_t total_read = 0;
  total_read += read(fd, options.msa_filename);
  total_read += read(fd, options.tree_filename);
  total_read += read(fd, options.prefix);
  total_read += read(fd, options.model_filename);
  total_read += read(fd, options.freqs_filename);
  total_read += read(fd, options.partition_filename);
  total_read += read(fd, options.data_type);
  total_read += read(fd, options.model_string);
  total_read += read(fd, options.rate_cats);
  total_read += read(fd, options.rate_category_types);
  total_read += read(fd, options.seed);
  total_read += read(fd, options.min_roots);
  total_read += read(fd, options.threads);
  total_read += read(fd, options.root_ratio);
  total_read += read(fd, options.abs_tolerance);
  total_read += read(fd, options.factor);
  total_read += read(fd, options.br_tolerance);
  total_read += read(fd, options.bfgs_tol);
  // total_read += read(fd, options.states);
  total_read += read(fd, options.silent);
  total_read += read(fd, options.exhaustive);
  total_read += read(fd, options.echo);
  total_read += read(fd, options.invariant_sites);
  total_read += read(fd, options.early_stop);
  debug_print(EMIT_LEVEL_DEBUG, "Read %lu bytes for cli_options", total_read);
  return total_read;
}

template <> size_t write(int fd, const partition_parameters_t &pp) {
  size_t total_written = 0;
  total_written += write(fd, pp.subst_rates);
  total_written += write(fd, pp.freqs);
  total_written += write(fd, pp.gamma_alpha);
  total_written += write(fd, pp.gamma_weights);
  return total_written;
}

template <> size_t read(int fd, partition_parameters_t &pp) {
  size_t total_read = 0;
  total_read += read(fd, pp.subst_rates);
  total_read += read(fd, pp.freqs);
  total_read += read(fd, pp.gamma_alpha);
  total_read += read(fd, pp.gamma_weights);
  return total_read;
}

template <>
std::pair<uint32_t, uint32_t> compute_checksum_components(
    const partition_parameters_t &pp, uint32_t a, uint32_t b) {
  return compute_checksum_components(
      a, b, pp.subst_rates, pp.freqs, pp.gamma_alpha, pp.gamma_weights);
}

checkpoint_t::checkpoint_t(const std::string &prefix) {
  _checkpoint_filename = prefix + ".ckp";
  _existing_results    = (access(_checkpoint_filename.c_str(), F_OK) != -1);
  _file_descriptor =
      open(_checkpoint_filename.c_str(), O_RDWR | O_APPEND | O_CREAT, 0640);
  if (_file_descriptor == -1) {
    throw std::runtime_error("Failed to open the checkpoint file");
  }
}

void checkpoint_t::clean() {
  if (!_existing_results) { return; }
  std::string backup_filename = _checkpoint_filename + ".bak";
  if (__MPI_RANK__ == 0) {
    auto lock = write_lock<fcntl_lock_behavior::block>();
    lseek(_file_descriptor, 0, SEEK_SET);
    auto copy_fd = open(
        backup_filename.c_str(), O_RDWR | O_CREAT | O_APPEND | O_EXCL, 0640);
    if (copy_fd == -1) {
      throw std::runtime_error(
          "Failed to open the new checkpoint when cleaning the checkpoint");
    }
    cli_options_t options;
    read_with_success(_file_descriptor, options);
    write_with_success(copy_fd, options);
    auto progress = read_results();
    for (auto result : progress) {
      write_with_checksum(copy_fd, result.first);
      write_with_checksum(copy_fd, result.second);
    }
    close(copy_fd);
    rename(backup_filename.c_str(), _checkpoint_filename.c_str());
  }
}

void checkpoint_t::write(const rd_result_t &result) {
  debug_print(
      EMIT_LEVEL_MPI_DEBUG, "Writing result with root id: %lu", result.root_id);
  write_with_checksum(_file_descriptor, result);
}

void checkpoint_t::write(
    const std::vector<partition_parameters_t> &parameters) {
  write_with_checksum(_file_descriptor, parameters);
}

void checkpoint_t::write(
    const rd_result_t &                        result,
    const std::vector<partition_parameters_t> &parameters) {
  auto lock = write_lock<fcntl_lock_behavior::block>();
  write(result);
  write(parameters);
}

void checkpoint_t::save_options(const cli_options_t &options) {
  if (!_existing_results) {
    auto lock = write_lock<fcntl_lock_behavior::block>();
    write_with_success(_file_descriptor, options);
  }
}

void checkpoint_t::load_options(cli_options_t &options) {
  if (_existing_results) {
    debug_string(EMIT_LEVEL_WARNING,
                 "Loading options from the checkpoint file");
    int read_fd = open(_checkpoint_filename.c_str(), O_RDONLY);
    read_with_success(read_fd, options);
    close(read_fd);
  }
}

std::vector<rd_result_t> checkpoint_t::current_progress() {
  std::vector<rd_result_t> results;
  auto                     progress = read_results();
  results.reserve(progress.size());
  for (auto &pair : progress) { results.push_back(pair.first); }
  return results;
}

std::vector<size_t> checkpoint_t::completed_indicies() {
  auto                results = current_progress();
  std::vector<size_t> completed_idx(results.size());
  for (size_t i = 0; i < results.size(); ++i) {
    completed_idx[i] = results[i].root_id;
  }
  return completed_idx;
}

bool checkpoint_t::existing_checkpoint() const { return _existing_results; }

int checkpoint_t::get_inode() {
  struct stat statbuf;

  int ret = fstat(_file_descriptor, &statbuf);
  if (ret == -1) {
    throw std::runtime_error{
        "There was an error getting the INODE of the checkpoint"};
  }
  return statbuf.st_ino;
}

checkpoint_t::~checkpoint_t() { close(_file_descriptor); }

checkpoint_t::checkpoint_t(checkpoint_t &&other) {
  _file_descriptor = other._file_descriptor;
}

checkpoint_t &checkpoint_t::operator=(checkpoint_t &&other) {
  _file_descriptor = other._file_descriptor;
  return *this;
}

void checkpoint_t::reload() {
  close(_file_descriptor);
  _file_descriptor =
      open(_checkpoint_filename.c_str(), O_RDWR | O_APPEND | O_CREAT, 0640);
  if (_file_descriptor == -1) {

    throw std::runtime_error{"Failed to reload the checkpoint file"};
  }
}

std::vector<std::pair<rd_result_t, std::vector<partition_parameters_t>>>
checkpoint_t::read_results() {
  std::vector<std::pair<rd_result_t, std::vector<partition_parameters_t>>>
       results;
  auto lock       = write_lock<fcntl_lock_behavior::block>();
  auto current_fd = fcntl(_file_descriptor, F_DUPFD, 0);

  lseek(current_fd, 0, SEEK_SET);

  // read and discard the options header to seek to the start of the results
  {
    cli_options_t tmp_opts;
    read_with_success(current_fd, tmp_opts);
  }

  auto current_position = lseek(current_fd, 0, SEEK_CUR);
  auto end_position     = lseek(current_fd, 0, SEEK_END);
  lseek(current_fd, current_position, SEEK_SET);

  while (current_position < end_position) {
    try {
      rd_result_t rdr;
      size_t      bytes_read = read_with_checksum(current_fd, rdr);
      if (bytes_read == expected_read_size<rd_result_t>()) {
      } else if (bytes_read == 0) {
        break;
      } else {
        throw std::runtime_error(
            "Unexpected read size when loading results from the checkpoint");
      }
      std::vector<partition_parameters_t> params;
      read_with_checksum(current_fd, params);
      results.push_back(std::make_pair(rdr, params));
      current_position = lseek(current_fd, 0, SEEK_CUR);
    } catch (checkpoint_read_success_failure &e) {
      debug_string(
          EMIT_LEVEL_WARNING,
          "Checkpoint file is corrupted, we will resume with what we can");
      break;
    }
  }
  close(current_fd);
  return results;
}

bool checkpoint_t::needs_cleaning() {
  auto lock       = write_lock<fcntl_lock_behavior::block>();
  auto current_fd = fcntl(_file_descriptor, F_DUPFD, 0);

  lseek(current_fd, 0, SEEK_SET);

  {
    cli_options_t tmp_opts;
    read_with_success(current_fd, tmp_opts);
  }

  auto current_position = lseek(current_fd, 0, SEEK_CUR);
  auto end_position     = lseek(current_fd, 0, SEEK_END);
  lseek(current_fd, current_position, SEEK_SET);
  while (current_position < end_position) {
    try {
      rd_result_t rdr;
      size_t      bytes_read = read_with_checksum(current_fd, rdr);
      if (bytes_read == expected_read_size<rd_result_t>()) {
      } else if (bytes_read == 0) {
        break;
      } else {
        throw std::runtime_error(
            "Unexpected read size when loading results from the checkpoint");
      }
      std::vector<partition_parameters_t> params;
      read_with_checksum(current_fd, params);
      current_position = lseek(current_fd, 0, SEEK_CUR);
    } catch (checkpoint_read_success_failure &e) {
      debug_string(
          EMIT_LEVEL_WARNING,
          "Checkpoint file is corrupted, we will resume with what we can");
      close(current_fd);
      return true;
    }
  }
  close(current_fd);
  return false;
}
