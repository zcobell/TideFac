/*------------------------------GPL---------------------------------------//
// This file is part of TideFac.
//
// (c) 2020 Zachary Cobell
//
// TideFac is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// TideFac is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with TideFac.  If not, see <http://www.gnu.org/licenses/>.
//------------------------------------------------------------------------*/
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
    i += 86400;
  }
}

static void tide_grid(benchmark::State &state) {
  TideFac t(-90.0, 90.0, 0.25);
  t.addMajor8();
  size_t i = 0;
  while (state.KeepRunning()) {
    t.calculate(i);
    i += 86400;
  }
}

BENCHMARK(tide_calc);
BENCHMARK(tide_grid);

BENCHMARK_MAIN();
