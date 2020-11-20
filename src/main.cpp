#include <cstdio>
#include <numeric>
extern "C" {
#include <libpll/pll.h>
}
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#ifdef MPI_BUILD
#include <mpi.h>
#endif
#include <cstdlib>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "debug.h"
int __VERBOSE__       = EMIT_LEVEL_PROGRESS;
int __MPI_RANK__      = 0;
int __MPI_NUM_TASKS__ = 1;
#include "checkpoint.hpp"
#include "model.hpp"
#include "msa.hpp"
#include "tree.hpp"
#include "util.hpp"

#define STRING(s) #s
#define STRINGIFY(s) STRING(s)
#define GIT_REV_STRING STRINGIFY(GIT_REV)
#define GIT_COMMIT_STRING STRINGIFY(GIT_COMMIT)
#define BUILD_DATE_STRING STRINGIFY(BUILD_DATE)

static void print_version() {
  debug_print(EMIT_LEVEL_IMPORTANT, "Version: %s", GIT_REV_STRING);
  debug_print(EMIT_LEVEL_IMPORTANT, "Build Commit: %s", GIT_COMMIT_STRING);
  debug_print(EMIT_LEVEL_IMPORTANT, "Build Date: %s", BUILD_DATE_STRING);
}

static void print_run_header(
    const std::chrono::time_point<std::chrono::system_clock> &start_time,
    uint64_t                                                  seed,
    size_t                                                    threads,
    int                                                       argv,
    char **                                                   argc) {
  time_t st = std::chrono::system_clock::to_time_t(start_time);
  char   time_string[256];
  std::strftime(time_string, sizeof(time_string), "%F %T", std::localtime(&st));
  debug_string(EMIT_LEVEL_IMPORTANT, "Running Root Digger");
  print_version();
  debug_print(EMIT_LEVEL_IMPORTANT, "Started: %s", time_string);
  debug_print(EMIT_LEVEL_IMPORTANT, "Seed: %lu", seed);
  debug_print(EMIT_LEVEL_IMPORTANT, "Number of threads per proc: %lu", threads);
#ifdef MPI_VERSION
  debug_print(EMIT_LEVEL_IMPORTANT, "Number of procs %d", __MPI_NUM_TASKS__);
#endif
  debug_print(EMIT_LEVEL_IMPORTANT,
              "Command: %s",
              combine_argv_argc(argv, argc).c_str());
  debug_string(EMIT_LEVEL_IMPORTANT,
               "Please report any bugs to "
               "https://groups.google.com/forum/#!forum/raxml");
}

static void print_usage() {
  if (__MPI_RANK__ != 0) { return; }
  std::cout
      << "Usage: rd [options]\n"
      << "Version: " << GIT_REV_STRING << "\n"
      << "Application Options:\n"
      << "  --msa [FILE]\n"
      << "         File containing the alignment.\n"
      << "  --tree [FILE]\n"
      << "         File containing the tree, with branch lengths.\n"
      << "  --partition [FILE]\n"
      << "         Optional file containing the partition specification.\n"
      << "         Format is the same as RAxML-NG partition file.\n"
      << "  --prefix [STRING]\n"
      << "         Prefix for the output files.\n"
      << "  --exhaustive\n"
      << "         Enable exhaustive mode. This will attempt to root a tree\n"
      << "         at every branch, and then report the results using LWR.\n"
      << "  --early-stop\n"
      << "         Enable early stopping. This will cause cause the search\n"
      << "         to terminate when the root placement is sufficently\n"
      << "         close for 2 consecutive iterations. How close they need\n"
      << "         to be can be controled by brtol. Is enabled by default for\n"
      << "         search mode and disabled by default for exhaustive mode.\n"
      << "  --no-early-stop\n"
      << "         Force disable early stop.\n"
      << "  --seed [NUMBER]\n"
      << "         Random seed to use. Optional\n"
      << "  --rate-cats [NUMBER]\n"
      << "         Number of rate categories to use for the model. Default\n"
      << "         is 1.\n"
      << "  --invariant-sites\n"
      << "         Enable invariant sites. Default is off.\n"
      << "  --min-roots [NUMBER]\n"
      << "         Minimum number of roots to start from. Optional,\n"
      << "         Default is 1.\n"
      << "  --root-ratio [NUMBER]\n"
      << "         Proportion of potential starting roots to attempt\n"
      << "         Default is 0.01\n"
      << "  --atol [NUMBER]\n"
      << "         Root optmization stopping tolerance. Increase this to \n"
      << "         improve results.Default is 1e-4\n"
      << "  --brtol [NUMBER]\n"
      << "         When early stop mode is enabled, this controls the\n"
      << "         distance required to trigger. Default is 1e-12\n"
      << "  --bfgstol [NUMBER]\n"
      << "         Tolerance for the BFGS steps. Default is 1e-7\n"
      << "  --factor [NUMBER]\n"
      << "         Factor for the BFGS steps. Default is 1e2\n"
      << "  --initial-root-strategy {random, midpoint, modified-mad}\n"
      << "         The strategy to pick the initial branches for rooting.\n"
      << "         Random is the default, and simply picks the branches at\n"
      << "         random. Midpoint uses the midpoint and similar branches to\n"
      << "         start. Modified MAD will use a modified version of mad to\n"
      << "         pick the starting branches. This can actually drive the\n"
      << "         so care should be taken when selecting this option.\n"
      << "         Default is random\n"
      << "  --threads [NUMBER]\n"
      << "         Number of threads to use\n"
      << "  --silent\n"
      << "         Suppress output except for the final tree\n"
      << "  --verbose\n"
      << "         Increase the verbosity level. Can be repeated to\n"
      << "         level further.\n"
      << "  --clean\n"
      << "         Clean the checkpoint file and exit. Normally, this should\n"
      << "         not be needed, but occasionally cleaining on a multi-node\n"
      << "         system can take a lot of time. In that case, use this flag\n"
      << "         on a single node, which will make RootDigger clean the\n"
      << "         checkpoint file so that the job can be run quickly on\n"
      << "         multi-node systems.\n"
      << std::endl;
}

cli_options_t parse_options(int argv, char **argc) {
  static struct option long_opts[] = {
      {"msa", required_argument, 0, 0},                   /* 0 */
      {"tree", required_argument, 0, 0},                  /* 1 */
      {"model", required_argument, 0, 0},                 /* 2 */
      {"seed", required_argument, 0, 0},                  /* 3 */
      {"verbose", no_argument, 0, 0},                     /* 4 */
      {"silent", no_argument, 0, 0},                      /* 5 */
      {"min-roots", required_argument, 0, 0},             /* 6 */
      {"root-ratio", required_argument, 0, 0},            /* 7 */
      {"atol", required_argument, 0, 0},                  /* 8 */
      {"brtol", required_argument, 0, 0},                 /* 9 */
      {"bfgstol", required_argument, 0, 0},               /* 10 */
      {"factor", required_argument, 0, 0},                /* 11 */
      {"partition", required_argument, 0, 0},             /* 12 */
      {"prefix", required_argument, 0, 0},                /* 13 */
      {"exhaustive", no_argument, 0, 0},                  /* 14 */
      {"early-stop", no_argument, 0, 0},                  /* 15 */
      {"no-early-stop", no_argument, 0, 0},               /* 16 */
      {"rate-cats", required_argument, 0, 0},             /* 17 */
      {"rate-cats-type", required_argument, 0, 0},        /* 18 */
      {"invariant-sites", no_argument, 0, 0},             /* 19 */
      {"initial-root-strategy", required_argument, 0, 0}, /* 20 */
      {"threads", required_argument, 0, 0},               /* 21 */
      {"version", no_argument, 0, 0},                     /* 22 */
      {"debug", no_argument, 0, 0},                       /* 23 */
      {"mpi-debug", no_argument, 0, 0},                   /* 24 */
      {"clean", no_argument, 0, 0},                       /* 25 */
      {"echo", no_argument, 0, 0},                        /* 26 */
      {"help", no_argument, 0, 0},                        /* 27 */
      {0, 0, 0, 0},
  };

  int c;
  int index = 0;
  optind    = 0;
  cli_options_t cli_options;
  while ((c = getopt_long_only(argv, argc, "", long_opts, &index)) == 0) {
    debug_print(EMIT_LEVEL_DEBUG, "parsing option index: %d", index);
    switch (index) {
    case 0: // msa
      cli_options.msa_filename = optarg;
      break;
    case 1: // tree
      cli_options.tree_filename = optarg;
      break;
    case 2: // model
      cli_options.model_string = optarg;
      break;
    case 3: // seed
      cli_options.seed = static_cast<uint64_t>(atol(optarg));
      break;
    case 4: // verbose
      __VERBOSE__ += 1;
      break;
    case 5: // silent
      __VERBOSE__        = 0;
      cli_options.silent = true;
      break;
    case 6: // min-roots
      cli_options.min_roots = static_cast<size_t>(atol(optarg));
      break;
    case 7: // root-ratio
      cli_options.root_ratio = atof(optarg);
      break;
    case 8: // atol
      cli_options.abs_tolerance = atof(optarg);
      break;
    case 9: // brtol
      cli_options.br_tolerance = atof(optarg);
      break;
    case 10: // bfgs_tol
      cli_options.bfgs_tol = atof(optarg);
      break;
    case 11: // factor
      cli_options.factor = atof(optarg);
      break;
    case 12: // partition
      cli_options.partition_filename = optarg;
      break;
    case 13: // prefix
      cli_options.prefix = optarg;
      break;
    case 14: // exhaustive
      cli_options.exhaustive = true;
      break;
    case 15: // early-stop
      cli_options.early_stop = initialized_flag_t::initialized_true;
      break;
    case 16: // no-early-stop
      cli_options.early_stop = initialized_flag_t::initialized_false;
      break;
    case 17: // rate-cats
      cli_options.rate_cats[0] = {(size_t)atol(optarg)};
      break;
    case 18: // rate-cats-type
      if (strcmp(optarg, "mean") == 0) {
        cli_options.rate_cats[0].rate_category_type = rate_category::MEAN;
      } else if (strcmp(optarg, "median") == 0) {
        cli_options.rate_cats[0].rate_category_type = rate_category::MEDIAN;
      } else if (strcmp(optarg, "free") == 0) {
        cli_options.rate_cats[0].rate_category_type = rate_category::FREE;
      }
      break;
    case 19: // invariant-sites
      cli_options.invariant_sites = true;
      break;
    case 20: // initial-root-strategy
      if (strcmp(optarg, "random") == 0) {
        cli_options.initial_root_strategy = {initial_root_strategy_t::random};
      } else if (strcmp(optarg, "midpoint") == 0) {
        cli_options.initial_root_strategy = {initial_root_strategy_t::midpoint};
      } else if (strcmp(optarg, "modified-mad") == 0) {
        cli_options.initial_root_strategy = {
            initial_root_strategy_t::modified_mad};
      } else {
        throw std::runtime_error{
            "An argument is required for --initial-root-strategy"};
      }

      break;
    case 21: // threads
      cli_options.threads = {(size_t)atol(optarg)};
      break;
    case 22: // version
      print_version();
      std::exit(0);
    case 23: // debug
      __VERBOSE__ = EMIT_LEVEL_DEBUG;
      break;
    case 24: // mpi-debug
      __VERBOSE__ = EMIT_LEVEL_MPI_DEBUG;
      break;
    case 25: // clean
      cli_options.clean = true;
      break;
    case 26: // echo
      cli_options.echo = true;
      break;
    case 27: // help
      print_usage();
      std::exit(0);
      break;
    case '?':
    case ':':
      print_usage();
      std::exit(0);
    default:
      throw std::invalid_argument("An argument was not recognized");
    }
  }
  if (c == '?' || c == ':') {
    print_usage();
    std::exit(1);
  }
#ifdef MPI_VERSION
  debug_string(EMIT_LEVEL_DEBUG, "Broadcasting seed");
  MPI_Bcast(&cli_options.seed, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
  debug_print(EMIT_LEVEL_MPI_DEBUG, "seed: %lu", cli_options.seed);
#endif

  if (cli_options.root_ratio < 0) {
    throw std::runtime_error("Root ratio is negative");
  }

  return cli_options;
}

/* This function defines what gets overridden by the checkpoint, and what
 * doesn't. There are plenty of things that need to stay the same, like the
 * search options, and the "model". But the number of threads/processes can
 * change. So this function handles that
 */
void merge_options_checkpoint(cli_options_t &cli_options,
                              checkpoint_t & checkpoint) {
  if (!checkpoint.existing_checkpoint()) { return; }

  cli_options_t checkpoint_options;
  checkpoint.load_options(checkpoint_options);
  checkpoint_options.threads = cli_options.threads;
  checkpoint_options.silent  = cli_options.silent;
  checkpoint_options.clean   = cli_options.clean;

  std::swap(cli_options, checkpoint_options);
}

void verify_options(const cli_options_t &cli_options) {
  if (cli_options.msa_filename.empty()) {
    std::cout << "No MSA was given, please supply an MSA" << std::endl;
    print_usage();
    std::exit(1);
  }
  if (cli_options.tree_filename.empty()) {
    std::cout << "No tree was given, please supply an tree" << std::endl;
    print_usage();
    std::exit(1);
  }
}

int wrapped_main(int argv, char **argc) {
#ifdef MPI_VERSION
  MPI_Comm_size(MPI_COMM_WORLD, &__MPI_NUM_TASKS__);
  MPI_Comm_rank(MPI_COMM_WORLD, &__MPI_RANK__);
#endif

#ifdef MPI_VERSION
  if (__MPI_NUM_TASKS__ == 1) {
    debug_string(EMIT_LEVEL_WARNING,
                 "Running MPI version with only 1 process, "
  }
#endif

  if (argv == 1) {
    print_usage();
    return 0;
  }

  auto start_time = std::chrono::system_clock::now();
  if (argv == 1) {
    print_usage();
    return 0;
  }
  try {

    auto cli_options = parse_options(argv, argc);
#ifdef MPI_VERSION
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    /* Use the tree path for the prefix */
    if (cli_options.prefix.empty()) {
      cli_options.prefix = cli_options.tree_filename;
    }

    checkpoint_t checkpoint(cli_options.prefix);
    merge_options_checkpoint(cli_options, checkpoint);
    if (__MPI_RANK__ == 0) {
      if (cli_options.clean) {
        debug_print(EMIT_LEVEL_IMPORTANT,
                    "Cleaning the checkpoint file %s",
                    checkpoint.get_filename().c_str());
        checkpoint.clean();
#ifdef MPI_VERSION
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        return 0;
      }
      checkpoint.save_options(cli_options);
      if (checkpoint.needs_cleaning()) { checkpoint.clean(); }
    } else if (cli_options.clean) {
#ifdef MPI_VERSION
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      return 0;
    }
#ifdef MPI_VERSION
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    checkpoint.reload();
    debug_print(
        EMIT_LEVEL_MPI_DEBUG, "Checkpoint inode %d", checkpoint.get_inode());

    if (cli_options.threads == 0) {
#ifdef MPI_VERSION
      cli_options.threads = 1;
#else
      cli_options.threads = sysutil_get_cpu_cores();
#endif
    }

#ifdef _OPENMP
    omp_set_num_threads(cli_options.threads);
#endif

    if (!cli_options.silent)
      print_run_header(
          start_time, cli_options.seed, cli_options.threads, argv, argc);

    debug_print(
        EMIT_LEVEL_INFO, "abs_tolerance: %.08f", cli_options.abs_tolerance);
    if (cli_options.exhaustive
        && cli_options.early_stop.convert_with_default(
            !cli_options.exhaustive)) {}

    constexpr const pll_state_t *map = pll_map_nt;

    /* Parse the model */
    if (!cli_options.model_string.empty()) {
      auto mi = parse_model_info(cli_options.model_string);
      cli_options.rate_cats.clear();

      cli_options.rate_cats.push_back(mi.ratehet_opts);

      if (mi.ratehet_opts.alpha_init) {
        debug_string(EMIT_LEVEL_WARNING,
                     "Ignoring alpha in model string as it currently "
                     "is not suported");
      }
      auto subst_str{mi.subst_str};

      for (auto &ch : subst_str) { ch = std::tolower(ch); }
      if (subst_str != "unrest") {
        debug_print(EMIT_LEVEL_WARNING,
                    "Ignoring subst matrix %s for model from command line"
                    ". Currently only UNREST is supported",
                    mi.subst_str.c_str());
      }
    }

    /* Parse the MSA */
    std::vector<msa_t> msa;
    msa_partitions_t   part_infos;
    if (cli_options.partition_filename.empty()) {
      msa.emplace_back(cli_options.msa_filename, map, cli_options.states);
    } else {
      msa_t unparted_msa{
          cli_options.msa_filename, map, cli_options.states, false};
      part_infos = parse_partition_file(cli_options.partition_filename);
      msa        = unparted_msa.partition(part_infos);
    }

    /* Parse the partitions */
    if (part_infos.size() > 0) {
      if (cli_options.rate_cats.size() > 0) {
        debug_string(EMIT_LEVEL_WARNING,
                     "Using rate categories from the partition file "
                     "over the option passed on the command line");
      }
      cli_options.rate_cats.clear();
      for (auto &p : part_infos) {
        auto rate_cats = p.model.ratehet_opts;
        if (rate_cats.rate_cats == 0) { rate_cats.rate_cats = 1; }
        cli_options.rate_cats.push_back(rate_cats);
        if (p.model.ratehet_opts.alpha_init) {
          debug_print(EMIT_LEVEL_WARNING,
                      "Ignoring alpha in partition %s as it currently "
                      "is not suported",
                      p.partition_name.c_str());
        }
        auto subst_str{p.model.subst_str};

        for (auto &ch : subst_str) { ch = std::tolower(ch); }

        if (subst_str != "unrest") {
          debug_print(EMIT_LEVEL_WARNING,
                      "Ignoring subst matrix %s for partition "
                      "%s. Currently only UNREST is supported",
                      p.model.subst_str.c_str(),
                      p.partition_name.c_str());
        }
      }
    }

    if (cli_options.rate_cats.size() == 1
        && cli_options.rate_cats[0].rate_cats == 0) {
      throw std::runtime_error("Rate categories cannot be zero");
    }

    /* Check the msa for validity */
    for (auto &m : msa) { m.valid_data(); }

    /* Make the tree */
    rooted_tree_t tree{cli_options.tree_filename};

    if (cli_options.min_roots > tree.root_count()) {
      throw std::runtime_error(
          "Min roots is larger than the number of roots on the tree");
    }

    model_t model{
        tree,
        msa,
        cli_options.rate_cats,
        cli_options.invariant_sites,
        cli_options.seed,
        cli_options.early_stop.convert_with_default(!cli_options.exhaustive)};
    try {
      model.initialize_partitions(msa);
    } catch (const invalid_empirical_frequencies_exception &) {
      model.initialize_partitions_uniform_freqs(msa);
    }

    if (cli_options.echo) { std::cout << tree.newick() << std::endl; }

    model.initialize();
    root_location_t final_rl;
    double          final_lh = -std::numeric_limits<double>::infinity();
    std::string     final_tree_string;
    if (!cli_options.exhaustive) {
      model.assign_indicies_by_rank_search(
          cli_options.min_roots,
          cli_options.root_ratio,
          static_cast<size_t>(__MPI_RANK__),
          static_cast<size_t>(__MPI_NUM_TASKS__),
          cli_options.initial_root_strategy,
          checkpoint);
#ifdef MPI_VERSION
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      auto tmp = model.search(cli_options.min_roots,
                              cli_options.root_ratio,
                              cli_options.abs_tolerance,
                              cli_options.bfgs_tol,
                              cli_options.br_tolerance,
                              cli_options.factor,
                              checkpoint);
      if (__MPI_RANK__ == 0) {
        model.finalize();
        final_rl = tmp.first;
        final_lh = tmp.second;
        std::ofstream outfile{cli_options.prefix + ".rooted.tree"};
        final_tree_string = model.rooted_tree(final_rl).newick(false);
        outfile << final_tree_string;
      }
    } else {

      model.assign_indicies_by_rank_exhaustive(
          static_cast<size_t>(__MPI_RANK__),
          static_cast<size_t>(__MPI_NUM_TASKS__),
          checkpoint);

#ifdef MPI_VERSION
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      auto tmp = model.exhaustive_search(cli_options.abs_tolerance,
                                         cli_options.bfgs_tol,
                                         cli_options.br_tolerance,
                                         cli_options.factor,
                                         checkpoint);
      if (__MPI_RANK__ == 0) {
        model.finalize();
        final_rl          = tmp.first;
        final_lh          = tmp.second;
        final_tree_string = model.virtual_rooted_tree(final_rl).newick();
        {
          std::ofstream outfile{cli_options.prefix + ".lwr.tree"};
          outfile << final_tree_string;
        }
        {
          std::ofstream outfile{cli_options.prefix + ".rooted.tree"};
          auto          tmp_tree = model.rooted_tree(final_rl);
          outfile << tmp_tree.newick(false);
        }
      }
    }

    if (!cli_options.silent) {
      debug_print(EMIT_LEVEL_IMPORTANT, "Final LogLH: %.5f", final_lh);
    }

    if (__MPI_RANK__ == 0) { std::cout << final_tree_string << std::endl; }

    auto end_time = std::chrono::system_clock::now();

    std::chrono::duration<double> duration = end_time - start_time;

    if (!cli_options.silent && __MPI_RANK__ == 0)
      std::cout << "Inference took: " << duration.count() << "s" << std::endl;

    if (__MPI_RANK__ == 0) {}
  } catch (const std::exception &e) {
    if (__MPI_RANK__ == 0) {
      std::cout << "There was an error during processing:\n"
                << e.what() << std::endl;
    }
    return 1;
  }

  return 0;
}

int main(int argv, char **argc) {
#ifdef MPI_VERSION
  MPI_Init(&argv, &argc);
#endif
  wrapped_main(argv, argc);
#ifdef MPI_VERSION
  MPI_Finalize();
#endif
}
