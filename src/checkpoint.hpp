#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "debug.h"
#include "tree.hpp"
#include "util.hpp"
#include <cctype>
#include <cstddef>
#include <exception>
#include <fcntl.h>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <unistd.h>
#include <vector>

struct rd_result_t {
  size_t root_id;
  double lh;
  double alpha;
};

typedef uint32_t field_flags_t;

#define CHECKPOINT_WRITE_SUCCESS_FLAG 1 << 0

class checkpoint_write_failure : std::runtime_error {
  using std::runtime_error::runtime_error;
};

class checkpoint_read_failure : std::runtime_error {
  using std::runtime_error::runtime_error;
};

class checkpoint_read_success_failure : checkpoint_read_failure {
  using checkpoint_read_failure::checkpoint_read_failure;
};

/* Implements the Adler-32 checksum
 * Implementation mostly taken from wikipedia
 */
template <typename T> uint32_t compute_checksum(const T &val) {
  constexpr const uint32_t MOD_ADLER = 65521;
  const uint8_t *data = reinterpret_cast<const uint8_t *>(&val);
  constexpr size_t data_len = sizeof(T) / sizeof(uint8_t);
  uint32_t a = 1;
  uint32_t b = 0;

  for (size_t i = 0; i < data_len; ++i) {
    a = (a + data[i]) % MOD_ADLER;
    b = b + a % MOD_ADLER;
  }
  return (b << 16) | a;
}

template <typename T> size_t write(int fd, const T &val) {
  auto result = write(fd, &val, sizeof(T));
  if (result < static_cast<ssize_t>(sizeof(T))) {
    throw checkpoint_write_failure{"Failed to write all data to the file"};
  }
  return static_cast<size_t>(result);
}

template <typename T> size_t write(int fd, const std::vector<T> &vals) {
  size_t result = 0;
  result += write(fd, vals.size());
  for (auto &v : vals) {
    result += write(fd, v);
  }
  return result;
}

template <typename T> size_t write_with_success(int fd, const T &val) {
  auto result = write(fd, val);
  field_flags_t fft = CHECKPOINT_WRITE_SUCCESS_FLAG;
  result += write(fd, fft);
  return result;
}

template <typename T> size_t write_with_checksum(int fd, const T &val) {
  auto bytes_written = write(fd, val);
  auto checksum = compute_checksum(val);
  auto pos = lseek(fd, 0, SEEK_CUR);
  debug_print(EMIT_LEVEL_MPI_DEBUG, "writing with checksum: %x to position %lu",
              checksum, pos);
  bytes_written += write(fd, checksum);
  return bytes_written;
}

template <typename T> size_t read(int fd, T &val) {
  auto result = read(fd, &val, sizeof(val));
  if (result < 0) {
    throw checkpoint_read_failure{"Failed to read a value"};
  }
  return static_cast<size_t>(result);
}

template <typename T> size_t read(int fd, std::vector<T> &vals) {
  std::vector<T> tmp_vals;

  size_t total_read = 0;

  size_t vector_size = 0;
  total_read += read(fd, vector_size);
  tmp_vals.reserve(vector_size);
  for (size_t i = 0; i < vector_size; ++i) {
    T tmp_val;
    total_read += read(fd, tmp_val);
    tmp_vals.push_back(tmp_val);
  }

  std::swap(tmp_vals, vals);
  return total_read;
}

template <typename T> size_t read_with_success(int fd, T &val) {
  T tmp_val;
  auto read_bytes = read(fd, tmp_val);
  field_flags_t ff;
  read_bytes += read(fd, ff);

  if (!(ff & CHECKPOINT_WRITE_SUCCESS_FLAG)) {
    throw checkpoint_read_success_failure{
        "The current read was unsuccessful due to an unsuccessful write flag"};
  }
  std::swap(val, tmp_val);
  return read_bytes;
}

template <typename T> size_t read_with_checksum(int fd, T &val) {
  T tmp_val;
  auto bytes_read = read(fd, tmp_val);
  uint32_t checksum = compute_checksum(tmp_val);
  uint32_t written_checksum;
  bytes_read += read(fd, written_checksum);

  if (checksum != written_checksum) {
    debug_print(EMIT_LEVEL_DEBUG, "computed checksum: %x, written checksum %x",
                checksum, written_checksum);
    throw checkpoint_read_success_failure{
        "The current read was unsuccessful due to an unsuccessful write flag"};
  }
  std::swap(val, tmp_val);
  return bytes_read;
}

template <typename T> constexpr size_t expected_read_size() {
  return sizeof(T) + sizeof(field_flags_t);
}

template <> constexpr size_t expected_read_size<rd_result_t>() {
  return sizeof(rd_result_t) + sizeof(uint32_t);
}

namespace fcntl_lock_behavior {
enum fcntl_lock_block_t {
  block,
  noblock,
};
}

template <fcntl_lock_behavior::fcntl_lock_block_t W> class fcntl_lock_t {
public:
  fcntl_lock_t(int file_descriptor, int lock_type)
      : _file_descriptor{fcntl(file_descriptor, F_DUPFD, 0)}, _file_lock{} {
    _file_lock.l_type = lock_type;
    auto ret = fcntl(_file_descriptor,
                     (W == fcntl_lock_behavior::block) ? F_SETLKW : F_SETLK,
                     &_file_lock);
    if (ret == -1) {
      throw std::runtime_error("failed to obtain the lock");
    }
  }

  fcntl_lock_t(const fcntl_lock_t &l) = delete;

  fcntl_lock_t(fcntl_lock_t &&l)
      : _file_descriptor{l._file_descriptor}, _file_lock{l._file_lock} {
    l._file_descriptor = -1;
  }

  ~fcntl_lock_t() {
    if (_file_descriptor != -1) {
      _file_lock.l_type = F_UNLCK;
      fcntl(_file_descriptor, F_SETLK, _file_lock);
      fsync(_file_descriptor);
      close(_file_descriptor);
    }
  }

  fcntl_lock_t &operator=(const fcntl_lock_t &l) = delete;

  fcntl_lock_t &operator=(fcntl_lock_t &&l) {
    _file_descriptor = fcntl(l._file_descriptor, F_DUPFD, 0);
    _file_lock = l._file_lock;
    l._file_descriptor = -1;
    return *this;
  }

private:
  int _file_descriptor;
  flock _file_lock;
};

class checkpoint_t {
public:
  checkpoint_t(const std::string &prefix);
  ~checkpoint_t();

  checkpoint_t(checkpoint_t &&);
  checkpoint_t &operator=(checkpoint_t &&);

  checkpoint_t(const checkpoint_t &) = delete;
  checkpoint_t &operator=(const checkpoint_t &) = delete;

  void clean();
  template <typename T> void write(const T &val) {
    auto lock = write_lock<fcntl_lock_behavior::block>();
    write_with_success(_file_descriptor, val);
  }
  void write(const rd_result_t &);
  void save_options(const cli_options_t &);
  void load_options(cli_options_t &);
  void reload();
  int get_inode();
  bool existing_checkpoint() const;
  std::vector<rd_result_t> current_progress();
  std::vector<size_t> completed_indicies();

private:
  template <fcntl_lock_behavior::fcntl_lock_block_t W>
  fcntl_lock_t<W> write_lock() {
    return fcntl_lock_t<W>(_file_descriptor, F_WRLCK);
  }

  std::string _checkpoint_filename;
  int _file_descriptor;
  bool _existing_results;
};

#endif
