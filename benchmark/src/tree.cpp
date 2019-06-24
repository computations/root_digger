#include <benchmark/benchmark.h>
#include <data.hpp>
#include <tree.hpp>

static void BM_tree_constructor(benchmark::State &state) {
  std::string filename = data_files_dna[state.range(0)].second;
  for (auto _ : state) {
    rooted_tree_t tree{filename};
  }
}

BENCHMARK(BM_tree_constructor)->Arg(0)->Arg(1);

static void BM_tree_reroot(benchmark::State &state) {
  rooted_tree_t tree{data_files_dna[state.range(0)].second};
  auto rl = tree.root_location(state.range(1));
  for (auto _ : state) {
    tree.root_by(rl);
  }
}

BENCHMARK(BM_tree_reroot)
    ->Args({0, 0})
    ->Args({0, 2})
    ->Args({1, 0})
    ->Args({1, 20})
    ->Args({1, 120});

static void BM_tree_generate_operations(benchmark::State &state) {
  rooted_tree_t tree{data_files_dna[state.range(0)].second};
  auto rl = tree.root_location(state.range(1));
  for (auto _ : state) {
    benchmark::DoNotOptimize(tree.generate_operations(rl));
  }
}

BENCHMARK(BM_tree_generate_operations)
    ->Args({0, 0})
    ->Args({0, 2})
    ->Args({1, 0})
    ->Args({1, 20})
    ->Args({1, 120});

static void BM_tree_generate_derivative_operations(benchmark::State &state) {
  rooted_tree_t tree{data_files_dna[state.range(0)].second};
  auto rl = tree.root_location(state.range(1));
  for (auto _ : state) {
    benchmark::DoNotOptimize(tree.generate_derivative_operations(rl));
  }
}

BENCHMARK(BM_tree_generate_derivative_operations)
    ->Args({0, 0})
    ->Args({0, 2})
    ->Args({1, 0})
    ->Args({1, 20})
    ->Args({1, 120});
