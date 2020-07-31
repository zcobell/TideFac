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
#include <iomanip>
#include <limits>

#include "catch.hpp"
#include "tidefac.h"

TEST_CASE("Major 8 Tide Factors", "[major8]") {
  TideFac f;
  f.addMajor8();
  f.calculate(Date(2011, 11, 1, 0, 0, 0), 29.7046);

  //...Amplitude
  REQUIRE(f.amplitude(0) == Approx(0.24289000000000002));
  REQUIRE(f.amplitude(1) == Approx(0.11342));
  REQUIRE(f.amplitude(2) == Approx(0.046544999999999996));
  REQUIRE(f.amplitude(3) == Approx(0.0307625));
  REQUIRE(f.amplitude(4) == Approx(0.14204250000000002));
  REQUIRE(f.amplitude(5) == Approx(0.1008475));
  REQUIRE(f.amplitude(6) == Approx(0.0193135));
  REQUIRE(f.amplitude(7) == Approx(0.04708));

  //...Equilibrium Argument
  REQUIRE(f.equilibriumArgument(0) == Approx(241.82802987977877));
  REQUIRE(f.equilibriumArgument(1) == Approx(359.88858705269575));
  REQUIRE(f.equilibriumArgument(2) == Approx(166.2219210330479));
  REQUIRE(f.equilibriumArgument(3) == Approx(97.65123191473414));
  REQUIRE(f.equilibriumArgument(4) == Approx(318.8245122078793));
  REQUIRE(f.equilibriumArgument(5) == Approx(278.40076786446616));
  REQUIRE(f.equilibriumArgument(6) == Approx(202.89697390456098));
  REQUIRE(f.equilibriumArgument(7) == Approx(50.63607437499934));

  //...Earth Tide Reduction Factor
  REQUIRE(f.earthTideReductionFactor(0) == Approx(0.693));
  REQUIRE(f.earthTideReductionFactor(1) == Approx(0.693));
  REQUIRE(f.earthTideReductionFactor(2) == Approx(0.693));
  REQUIRE(f.earthTideReductionFactor(3) == Approx(0.693));
  REQUIRE(f.earthTideReductionFactor(4) == Approx(0.736));
  REQUIRE(f.earthTideReductionFactor(5) == Approx(0.695));
  REQUIRE(f.earthTideReductionFactor(6) == Approx(0.695));
  REQUIRE(f.earthTideReductionFactor(7) == Approx(0.706));

  //...Frequency
  REQUIRE(f.frequency(0) == Approx(0.00014049900478554353));
  REQUIRE(f.frequency(1) == Approx(0.00014538592669112763));
  REQUIRE(f.frequency(2) == Approx(0.00013788101090755203));
  REQUIRE(f.frequency(3) == Approx(0.00014590952546672593));
  REQUIRE(f.frequency(4) == Approx(7.295476273336297e-05));
  REQUIRE(f.frequency(5) == Approx(6.754424205218055e-05));
  REQUIRE(f.frequency(6) == Approx(6.492624817418905e-05));
  REQUIRE(f.frequency(7) == Approx(7.26056968829641e-05));

  //...Node Factor
  REQUIRE(f.nodeFactor(0) == Approx(1.0088910309243924));
  REQUIRE(f.nodeFactor(1) == Approx(0.9994881350641318));
  REQUIRE(f.nodeFactor(2) == Approx(1.010033667119077));
  REQUIRE(f.nodeFactor(3) == Approx(0.9493529816818805));
  REQUIRE(f.nodeFactor(4) == Approx(0.9865471373144313));
  REQUIRE(f.nodeFactor(5) == Approx(0.9761350226537572));
  REQUIRE(f.nodeFactor(6) == Approx(0.9748772611466947));
  REQUIRE(f.nodeFactor(7) == Approx(1.0016946305592558));

  //...Nodefactor Correction
  REQUIRE(f.nodefactorCorrection(0) == Approx(2.1731074244177875));
  REQUIRE(f.nodefactorCorrection(1) == Approx(359.88858705269575));
  REQUIRE(f.nodefactorCorrection(2) == Approx(1.8970867762275152));
  REQUIRE(f.nodefactorCorrection(3) == Approx(17.769348613397614));
  REQUIRE(f.nodefactorCorrection(4) == Approx(8.883570557211046));
  REQUIRE(f.nodefactorCorrection(5) == Approx(348.68678705977345));
  REQUIRE(f.nodefactorCorrection(6) == Approx(348.5130812984089));
  REQUIRE(f.nodefactorCorrection(7) == Approx(0.5770160256675876));

  //...Astronomical Argument
  REQUIRE(f.astronomicArgument(0) == Approx(239.65492245536097));
  REQUIRE(f.astronomicArgument(1) == Approx(0.0));
  REQUIRE(f.astronomicArgument(2) == Approx(164.32483425682037));
  REQUIRE(f.astronomicArgument(3) == Approx(79.88188330133653));
  REQUIRE(f.astronomicArgument(4) == Approx(309.94094165066826));
  REQUIRE(f.astronomicArgument(5) == Approx(289.7139808046927));
  REQUIRE(f.astronomicArgument(6) == Approx(214.3838926061521));
  REQUIRE(f.astronomicArgument(7) == Approx(50.059058349331735));
}

TEST_CASE("Major 8 Mean Values", "[major8mean]") {
  TideFac f;
  f.addMajor8();
  f.calculate(Date(2011, 11, 1, 0, 0, 0), Date(2011, 12, 1, 0, 0, 0), 29.7046);

  //...Amplitude
  REQUIRE(f.amplitude(0) == Approx(0.24289000000000002));
  REQUIRE(f.amplitude(1) == Approx(0.11342));
  REQUIRE(f.amplitude(2) == Approx(0.046544999999999996));
  REQUIRE(f.amplitude(3) == Approx(0.0307625));
  REQUIRE(f.amplitude(4) == Approx(0.14204250000000002));
  REQUIRE(f.amplitude(5) == Approx(0.1008475));
  REQUIRE(f.amplitude(6) == Approx(0.0193135));
  REQUIRE(f.amplitude(7) == Approx(0.04708));

  //...Equilibrium Argument
  REQUIRE(f.equilibriumArgument(0) == Approx(241.81892023876082));
  REQUIRE(f.equilibriumArgument(1) == Approx(359.88912769059374));
  REQUIRE(f.equilibriumArgument(2) == Approx(166.21154687086076));
  REQUIRE(f.equilibriumArgument(3) == Approx(97.60373664907745));
  REQUIRE(f.equilibriumArgument(4) == Approx(318.81304101465037));
  REQUIRE(f.equilibriumArgument(5) == Approx(278.3944623654454));
  REQUIRE(f.equilibriumArgument(6) == Approx(202.89461836204345));
  REQUIRE(f.equilibriumArgument(7) == Approx(50.632943102366085));

  //...Earth Tide Reduction Factor
  REQUIRE(f.earthTideReductionFactor(0) == Approx(0.693));
  REQUIRE(f.earthTideReductionFactor(1) == Approx(0.693));
  REQUIRE(f.earthTideReductionFactor(2) == Approx(0.693));
  REQUIRE(f.earthTideReductionFactor(3) == Approx(0.693));
  REQUIRE(f.earthTideReductionFactor(4) == Approx(0.736));
  REQUIRE(f.earthTideReductionFactor(5) == Approx(0.695));
  REQUIRE(f.earthTideReductionFactor(6) == Approx(0.695));
  REQUIRE(f.earthTideReductionFactor(7) == Approx(0.706));

  //...Frequency
  REQUIRE(f.frequency(0) == Approx(0.00014049900478554353));
  REQUIRE(f.frequency(1) == Approx(0.00014538592669112763));
  REQUIRE(f.frequency(2) == Approx(0.00013788101090755203));
  REQUIRE(f.frequency(3) == Approx(0.00014590952546672593));
  REQUIRE(f.frequency(4) == Approx(7.295476273336297e-05));
  REQUIRE(f.frequency(5) == Approx(6.754424205218055e-05));
  REQUIRE(f.frequency(6) == Approx(6.492624817418905e-05));
  REQUIRE(f.frequency(7) == Approx(7.26056968829641e-05));

  //...Node Factor
  REQUIRE(f.nodeFactor(0) == Approx(1.0093462359746201));
  REQUIRE(f.nodeFactor(1) == Approx(0.9994506447765682));
  REQUIRE(f.nodeFactor(2) == Approx(1.0102272459290909));
  REQUIRE(f.nodeFactor(3) == Approx(0.945622399462162));
  REQUIRE(f.nodeFactor(4) == Approx(0.9849031936868394));
  REQUIRE(f.nodeFactor(5) == Approx(0.9738428648398829));
  REQUIRE(f.nodeFactor(6) == Approx(0.9719888798690001));
  REQUIRE(f.nodeFactor(7) == Approx(1.0019341617997175));

  //...Nodefactor Correction
  REQUIRE(f.nodefactorCorrection(0) == Approx(2.163997783399843));
  REQUIRE(f.nodefactorCorrection(1) == Approx(359.88912769059374));
  REQUIRE(f.nodefactorCorrection(2) == Approx(1.8867126140403705));
  REQUIRE(f.nodefactorCorrection(3) == Approx(17.721853347740925));
  REQUIRE(f.nodefactorCorrection(4) == Approx(8.87209936398213));
  REQUIRE(f.nodefactorCorrection(5) == Approx(348.6804815607527));
  REQUIRE(f.nodefactorCorrection(6) == Approx(348.5107257558914));
  REQUIRE(f.nodefactorCorrection(7) == Approx(0.5738847530343619));

  //...Astronomical Argument
  REQUIRE(f.astronomicArgument(0) == Approx(239.65492245536097));
  REQUIRE(f.astronomicArgument(1) == Approx(0.0));
  REQUIRE(f.astronomicArgument(2) == Approx(164.32483425682037));
  REQUIRE(f.astronomicArgument(3) == Approx(79.88188330133653));
  REQUIRE(f.astronomicArgument(4) == Approx(309.94094165066826));
  REQUIRE(f.astronomicArgument(5) == Approx(289.7139808046927));
  REQUIRE(f.astronomicArgument(6) == Approx(214.3838926061521));
  REQUIRE(f.astronomicArgument(7) == Approx(50.059058349331735));
}
