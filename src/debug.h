#ifndef RD_DEBUG
#define RD_DEBUG

#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <sys/ioctl.h>
#include <time.h>
#include <unistd.h>

const auto CLOCK_START = std::chrono::high_resolution_clock::now();
extern int __VERBOSE__;
extern int __MPI_RANK__;
extern int __MPI_NUM_TASKS__;

#ifdef RD_DEBUG_FLAG
#define DEBUG_IF_FLAG 1
#else
#define DEBUG_IF_FLAG 0
#endif

#ifdef RD_DEBUG
#define RD_DEBUG_ASSERT_FLAG 1
#else
#define RD_DEBUG_ASSERT_FLAG 0
#endif

#define EMIT_LEVEL_IMPORTANT 0
#define EMIT_LEVEL_ERROR 1
#define EMIT_LEVEL_WARNING 2
#define EMIT_LEVEL_PROGRESS 3
#define EMIT_LEVEL_MPROGRESS 4
#define EMIT_LEVEL_INFO 5
#define EMIT_LEVEL_DEBUG 6
#define EMIT_LEVEL_MPI_DEBUG 7

#define progress_macro(i, k)                                                   \
  (((std::chrono::high_resolution_clock::now() - CLOCK_START).count() /        \
    static_cast<double>(i)) *                                                  \
   (static_cast<double>(k - i)) / 1e9 / 3600.0)

#define print_clock                                                            \
  do {                                                                         \
    std::chrono::duration<double> diff =                                       \
        std::chrono::high_resolution_clock::now() - CLOCK_START;               \
    printf("[%.2f] ", diff.count());                                           \
  } while (0)

#define debug_print(level, fmt, ...)                                           \
  do {                                                                         \
    if (DEBUG_IF_FLAG || level < EMIT_LEVEL_DEBUG) {                           \
      if (__VERBOSE__ >= level &&                                              \
          (__MPI_RANK__ == 0 || level == EMIT_LEVEL_MPI_DEBUG)) {              \
        print_clock;                                                           \
        if (__VERBOSE__ >= EMIT_LEVEL_DEBUG) {                                 \
          printf("[%s:%d]", __func__, __LINE__);                               \
          if (level == EMIT_LEVEL_MPI_DEBUG) {                                 \
            printf(" [Rank: %d]", __MPI_RANK__);                               \
          }                                                                    \
          printf(" : ");                                   \
        }                                                                      \
        if (level == EMIT_LEVEL_WARNING) {                                     \
          printf("[Warning] ");                                                \
        }                                                                      \
        if (level == EMIT_LEVEL_ERROR) {                                       \
          printf("[ERROR] ");                                                  \
        }                                                                      \
        printf(fmt "\n", __VA_ARGS__);                                         \
      }                                                                        \
    }                                                                          \
  } while (0)

#define debug_string(level, x)                                                 \
  do {                                                                         \
    debug_print(level, "%s", x);                                               \
  } while (0)

#define print_trace()                                                          \
  do {                                                                         \
    if (DEBUG_IF_FLAG) {                                                       \
      void *callstack[128];                                                    \
      int frames = backtrace(callstack, 128);                                  \
      char **bt_symbols = backtrace_symbols(callstack, frames);                \
      print_clock;                                                             \
      printf("BACKTRACE AT %s:%d:%s():\n", __FILE__, __LINE__,        \
              __func__);                                                       \
      for (int i = 0; i < frames; ++i) {                                       \
        print_clock;                                                           \
        printf("%s\n", bt_symbols[i]);                                \
      }                                                                        \
    }                                                                          \
  } while (0)

#define assert_string(cond, comment)                                           \
  do {                                                                         \
    if (DEBUG_IF_FLAG) {                                                       \
      {                                                                        \
        if (!(cond)) {                                                         \
          print_clock;                                                         \
          fprintf(                                                             \
              stderr,                                                          \
              "assertion \"%s\" failed: file: %s, line: %d, comment: %s\n",    \
              #cond, __FILE__, __LINE__, comment);                             \
          abort();                                                             \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  } while (0)
#endif
