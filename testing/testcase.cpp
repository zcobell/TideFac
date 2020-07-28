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
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "tidefac.h"

TEST_CASE("Major 8 Tide Factors", "[major8]") {
  TideFac f;
  f.addMajor8();
  f.calculate(Date(2011, 11, 1, 0, 0, 0), 29.7046);

  //...Amplitude
  REQUIRE(f.amplitude(0) == 0.24289000000000002);
  REQUIRE(f.amplitude(1) == 0.11342);
  REQUIRE(f.amplitude(2) == 0.046544999999999996);
  REQUIRE(f.amplitude(3) == 0.0307625);
  REQUIRE(f.amplitude(4) == 0.14204250000000002);
  REQUIRE(f.amplitude(5) == 0.1008475);
  REQUIRE(f.amplitude(6) == 0.0193135);
  REQUIRE(f.amplitude(7) == 0.04708);

  //...Equilibrium Argument
  REQUIRE(f.equilibriumArgument(0) == 241.82802987977877);
  REQUIRE(f.equilibriumArgument(1) == 359.88858705269575);
  REQUIRE(f.equilibriumArgument(2) == 166.2219210330479);
  REQUIRE(f.equilibriumArgument(3) == 97.65123191473414);
  REQUIRE(f.equilibriumArgument(4) == 318.8245122078793);
  REQUIRE(f.equilibriumArgument(5) == 278.40076786446616);
  REQUIRE(f.equilibriumArgument(6) == 202.89697390456098);
  REQUIRE(f.equilibriumArgument(7) == 50.63607437499934);

  //...Earth Tide Reduction Factor
  REQUIRE(f.earthTideReductionFactor(0) == 0.693);
  REQUIRE(f.earthTideReductionFactor(1) == 0.693);
  REQUIRE(f.earthTideReductionFactor(2) == 0.693);
  REQUIRE(f.earthTideReductionFactor(3) == 0.693);
  REQUIRE(f.earthTideReductionFactor(4) == 0.736);
  REQUIRE(f.earthTideReductionFactor(5) == 0.695);
  REQUIRE(f.earthTideReductionFactor(6) == 0.695);
  REQUIRE(f.earthTideReductionFactor(7) == 0.706);

  //...Frequency
  REQUIRE(f.frequency(0) == 0.00014049900478554353);
  REQUIRE(f.frequency(1) == 0.00014538592669112763);
  REQUIRE(f.frequency(2) == 0.00013788101090755203);
  REQUIRE(f.frequency(3) == 0.00014590952546672593);
  REQUIRE(f.frequency(4) == 7.295476273336297e-05);
  REQUIRE(f.frequency(5) == 6.754424205218055e-05);
  REQUIRE(f.frequency(6) == 6.492624817418905e-05);
  REQUIRE(f.frequency(7) == 7.26056968829641e-05);

  //...Node Factor
  REQUIRE(f.nodeFactor(0) == 1.0088910309243924);
  REQUIRE(f.nodeFactor(1) == 0.9994881350641318);
  REQUIRE(f.nodeFactor(2) == 1.010033667119077);
  REQUIRE(f.nodeFactor(3) == 0.9493529816818805);
  REQUIRE(f.nodeFactor(4) == 0.9865471373144313);
  REQUIRE(f.nodeFactor(5) == 0.9761350226537572);
  REQUIRE(f.nodeFactor(6) == 0.9748772611466947);
  REQUIRE(f.nodeFactor(7) == 1.0016946305592558);


}


