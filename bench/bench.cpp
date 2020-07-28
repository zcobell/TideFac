
#include "benchmark/benchmark.h"
#include "tidefac.h"

static void tide_calc(benchmark::State &state) {
  TideFac t;
  Date refDate(2011, 1, 1, 0, 0, 0);
  t.addMajor8();
  t.setRefTime(refDate);
  size_t i = 0;
  while (state.KeepRunning()) {
    t.calculate(i, 29.5);
    i = i * 86400;
  }
}

BENCHMARK(tide_calc);

BENCHMARK_MAIN();
