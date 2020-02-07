#include <benchmark/benchmark.h>
#include <data.hpp>
#include <model.hpp>

static void BM_model_constructor(benchmark::State &state) {
  std::vector<msa_t> msa;
  size_t data_index = static_cast<size_t>(state.range(0));
  msa.emplace_back(data_files_dna[data_keys[data_index]].first);
  rooted_tree_t tree{data_files_dna[data_keys[data_index]].second};
  uint32_t seed = (uint32_t)std::rand();
  model_params_t freqs{.25, .25, .25, .25};
  for (auto _ : state) {
    model_t model{tree, msa, {1}, true, seed, false};
  }
}

BENCHMARK(BM_model_constructor)->Arg(0lu)->Arg(1lu);

static void BM_LH_computation(benchmark::State &state) {
  std::vector<msa_t> msa;
  size_t data_index = static_cast<size_t>(state.range(0));
  msa.emplace_back(data_files_dna[data_keys[data_index]].first);
  rooted_tree_t tree{data_files_dna[data_keys[data_index]].second};
  uint32_t seed = (uint32_t)std::rand();
  model_params_t freqs{.25, .25, .25, .25};
  model_t model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  auto rl = tree.root_location(static_cast<size_t>(state.range(1)));
  for (auto _ : state) {
    benchmark::DoNotOptimize(model.compute_lh(rl));
  }
}

BENCHMARK(BM_LH_computation)
    ->Args({0lu, 0})
    ->Args({0lu, 2})
    ->Args({1lu, 0})
    ->Args({1lu, 20})
    ->Args({1lu, 120});

static void BM_DLH_computation(benchmark::State &state) {
  std::vector<msa_t> msa;
  size_t data_index = static_cast<size_t>(state.range(0));
  msa.emplace_back(data_files_dna[data_keys[data_index]].first);
  rooted_tree_t tree{data_files_dna[data_keys[data_index]].second};
  uint32_t seed = (uint32_t)std::rand();
  model_params_t freqs{.25, .25, .25, .25};
  model_t model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  auto rl = tree.root_location(static_cast<size_t>(state.range(1)));
  model.compute_lh(rl);
  for (auto _ : state) {
    benchmark::DoNotOptimize(model.compute_dlh(rl));
  }
}

BENCHMARK(BM_DLH_computation)
    ->Args({0lu, 0})
    ->Args({0lu, 2})
    ->Args({1lu, 0})
    ->Args({1lu, 20})
    ->Args({1lu, 120});

static void BM_LH_root_computation(benchmark::State &state) {
  std::vector<msa_t> msa;
  msa.emplace_back(data_files_dna["101.phy"].first);
  rooted_tree_t tree{data_files_dna["101.phy"].second};
  uint32_t seed = (uint32_t)std::rand();
  model_t model{tree, msa, {1}, true, seed, false};
  model.initialize_partitions_uniform_freqs(msa);
  auto rl = tree.root_location(static_cast<size_t>(state.range(0)));
  model.compute_lh(rl);
  for (auto _ : state) {
    benchmark::DoNotOptimize(model.compute_lh_root(rl));
  }
}

BENCHMARK(BM_LH_root_computation)->Arg(0)->Arg(20)->Arg(120);
