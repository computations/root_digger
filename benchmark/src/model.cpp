#include <benchmark/benchmark.h>
#include <data.hpp>
#include <model.hpp>

static void BM_model_constructor(benchmark::State &state) {
  msa_t msa{data_files_dna[state.range(0)].first};
  rooted_tree_t tree{data_files_dna[state.range(0)].second};
  uint64_t seed = std::rand();
  model_params_t freqs{.25, .25, .25, .25};
  for (auto _ : state) {
    model_t{tree, msa, freqs, seed};
  }
}

BENCHMARK(BM_model_constructor)->Arg(0)->Arg(1);

static void BM_LH_computation(benchmark::State &state) {
  msa_t msa{data_files_dna[state.range(0)].first};
  rooted_tree_t tree{data_files_dna[state.range(0)].second};
  uint64_t seed = std::rand();
  model_params_t freqs{.25, .25, .25, .25};
  model_t model{tree, msa, freqs, seed};
  auto rl = tree.root_location(state.range(1));
  for (auto _ : state) {
    benchmark::DoNotOptimize(model.compute_lh(rl));
  }
}

BENCHMARK(BM_LH_computation)
    ->Args({0, 0})
    ->Args({0, 2})
    ->Args({1, 0})
    ->Args({1, 20})
    ->Args({1, 120});

static void BM_DLH_computation(benchmark::State &state) {
  msa_t msa{data_files_dna[state.range(0)].first};
  rooted_tree_t tree{data_files_dna[state.range(0)].second};
  uint64_t seed = std::rand();
  model_params_t freqs{.25, .25, .25, .25};
  model_t model{tree, msa, freqs, seed};
  auto rl = tree.root_location(state.range(1));
  for (auto _ : state) {
    benchmark::DoNotOptimize(model.compute_dlh(rl));
  }
}

BENCHMARK(BM_DLH_computation)
    ->Args({0, 0})
    ->Args({0, 2})
    ->Args({1, 0})
    ->Args({1, 20})
    ->Args({1, 120});

static void BM_LH_root_computation(benchmark::State &state) {
  msa_t msa{data_files_dna[1].first};
  rooted_tree_t tree{data_files_dna[1].second};
  uint64_t seed = std::rand();
  model_t model{tree, msa, seed};
  auto rl = tree.root_location(state.range(0));
  model.compute_lh(rl);
  for (auto _ : state) {
    benchmark::DoNotOptimize(model.compute_lh_root(rl));
  }
}

BENCHMARK(BM_LH_root_computation)->Arg(0)->Arg(20)->Arg(120);
