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
#ifndef TIDEFAC_H
#define TIDEFAC_H

#include <array>
#include <complex>
#include <vector>

#include "TideFac_global.h"
#include "date.h"

class TideFac {
 public:
  TIDEFAC_EXPORT TideFac();
  TIDEFAC_EXPORT TideFac(const Date &refTime);

  int TIDEFAC_EXPORT addConstituent(const std::string &harmonic);
  int TIDEFAC_EXPORT addConstituent(const char *harmonic);
  int TIDEFAC_EXPORT addMajor8();

  void TIDEFAC_EXPORT setRefTime(const Date &refTime);
  Date TIDEFAC_EXPORT refTime() const;

  void TIDEFAC_EXPORT calculate(const Date &d, const double latitude);
  void TIDEFAC_EXPORT calculate(const size_t dt, const double latitude);
  void TIDEFAC_EXPORT calculate(const double dt, const double latitude);

  void TIDEFAC_EXPORT show() const;

  Date TIDEFAC_EXPORT curTime() const;

  std::string TIDEFAC_EXPORT name(size_t index);
  double TIDEFAC_EXPORT amplitude(size_t index);
  double TIDEFAC_EXPORT frequency(size_t index);
  double TIDEFAC_EXPORT earthTideReductionFactor(size_t index);
  double TIDEFAC_EXPORT nodeFactor(size_t index);
  double TIDEFAC_EXPORT equilibriumArgument(size_t index);

 private:
  struct Tide {
    std::string name;
    double amp;
    double freq;
    double etrf;
    double nodefactor;
    double eqarg;
  };

  static constexpr double pi() { return 3.14159265358979323846264338327950288; }
  static constexpr double twopi() { return pi() * 2.0; }
  static constexpr double pi180() { return pi() / 180.0; }

  std::tuple<std::vector<double>, std::vector<double>, std::vector<double> >
  computeAstronomicalArguments(const Date &d, const double latitude);

  void computeOrbitalParameters(const Date &d, std::array<double, 6> *astro,
                                std::array<double, 6> *ader = nullptr);

  Date m_refTime;
  Date m_curTime;
  std::vector<std::string> m_constituentNames;
  std::vector<size_t> m_constituentIndex;
  std::vector<Tide> m_tides;
};

#endif  // TIDEFAC_H
