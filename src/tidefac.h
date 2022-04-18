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
#include <string>
#include <utility>
#include <vector>

#include "date.h"
#include "tidefac_global.h"

/**
 * @class TideFac
 * @author Zachary Cobell
 * @brief Class to handle the calculation of tide parameters on the earth
 * @copyright Copyright 2015-2020 Zachary Cobell. All Rights Reserved. This
 * project is released under the terms of the GNU General Public License v3
 *
 * The TideFac class is designed to compute the tide factors for use in
 * numerical modeling. The methodology is developed from UTide (Codiga, 2011)
 *
 */
class TideFac {
 public:
  TIDEFAC_EXPORT TideFac();
  TIDEFAC_EXPORT explicit TideFac(const Date &refTime);
  TIDEFAC_EXPORT TideFac(double latmin, double latmax, double resolution);

  int TIDEFAC_EXPORT generateLatitudeGrid(double latmin, double latmax,
                                          double resolution);

  int TIDEFAC_EXPORT addConstituent(std::string harmonic);
  int TIDEFAC_EXPORT addConstituent(const char *harmonic);
  int TIDEFAC_EXPORT addMajor8();

  void TIDEFAC_EXPORT setRefTime(const Date &refTime);
  Date TIDEFAC_EXPORT refTime() const;

  void TIDEFAC_EXPORT calculate(const Date &d, double latitude);
  void TIDEFAC_EXPORT calculate(size_t dt, double latitude);
  void TIDEFAC_EXPORT calculate(double dt, double latitude);
  void TIDEFAC_EXPORT calculate(const Date &d1, const Date &d2,
                                double latitude);
  void TIDEFAC_EXPORT calculate(const Date &d);
  void TIDEFAC_EXPORT calculate(size_t dt);
  void TIDEFAC_EXPORT calculate(double dt);

  void TIDEFAC_EXPORT show(size_t index = 0) const;

  Date TIDEFAC_EXPORT curTime() const;

  std::string TIDEFAC_EXPORT name(size_t index, size_t gridIndex = 0) const;
  double TIDEFAC_EXPORT amplitude(size_t index, size_t gridIndex = 0) const;
  double TIDEFAC_EXPORT frequency(size_t index, size_t gridIndex = 0) const;
  double TIDEFAC_EXPORT earthTideReductionFactor(size_t index,
                                                 size_t gridIndex = 0) const;
  double TIDEFAC_EXPORT nodeFactor(size_t index, size_t gridIndex = 0) const;
  double TIDEFAC_EXPORT equilibriumArgument(size_t index, size_t gridIndex = 0);
  double TIDEFAC_EXPORT nodefactorCorrection(size_t index,
                                             size_t gridIndex = 0) const;
  double TIDEFAC_EXPORT astronomicArgument(size_t index,
                                           size_t gridIndex = 0) const;

  int TIDEFAC_EXPORT getInterpolationFactors(double latitude,
                                             unsigned long &gridIndex,
                                             double &weight) const;

 private:
  struct Tide {
    Tide(std::string name, const double amp, const double freq,
         const double etrf, const double nodefactor, const double eqarg,
         const double nodecorrection, const double astroarg)
        : name(std::move(name)),
          amp(amp),
          freq(freq),
          etrf(etrf),
          nodefactor(nodefactor),
          eqarg(eqarg),
          nodecorrection(nodecorrection),
          astroarg(astroarg) {}
    std::string name;
    double amp;
    double freq;
    double etrf;
    double nodefactor;
    double eqarg;
    double nodecorrection;
    double astroarg;
  };

  struct OrbitalDateData {
    long seconds;
    Date julian_date;
    double days;
    double dd;
    double nd;
    std::array<double, 4> args;
  };

  static void toUpper(std::string &s);

  static constexpr double pi() { return M_PI; }
  static constexpr double twopi() { return pi() * 2.0; }
  static constexpr double deg2rad() { return pi() / 180.0; }
  static constexpr double rad2deg() { return 180.0 / pi(); }

  template <typename T>
  static constexpr T meanNumTides() {
    return static_cast<T>(1000.0);
  }

  static std::tuple<std::vector<double>, std::vector<double>,
                    std::vector<double>>
  computeAstronomicalArguments(const std::array<double, 6> &astro,
                               double latitude);

  static std::array<double, 162> computeUU(const std::array<double, 6> &astro);

  static std::array<std::complex<double>, 162> computeMat(
      const std::array<double, 162> &uu, const std::array<double, 162> &rr);

  static std::array<double, 162> computeRR(double latitude);

  OrbitalDateData computeOrbitalDates(const Date &d);

  std::array<double, 6> computeSomeOrbitalParameters(const Date &d);
  static std::array<double, 6> computeSomeOrbitalParameters(const OrbitalDateData &d);

  std::tuple<std::array<double, 6>, std::array<double, 6>>
  computeAllOrbitalParameters(const Date &d);

  static void computePrimaryFactors(
      std::vector<std::complex<double>> &F, std::vector<double> &U,
      std::vector<double> &V, const std::array<std::complex<double>, 162> &mat,
      const std::array<double, 6> &astro);

  static void computeMinorCorrections(std::vector<std::complex<double>> &F,
                                      std::vector<double> &U,
                                      std::vector<double> &V);

  static std::vector<double> complexToReal(
      std::vector<std::complex<double>> &complex);

  static std::vector<double> generateGrid(double min, double max, double res);

  static void reorientAngularParameters(double &nc, double &eq, double &aa);

  TideFac::Tide generateTide(size_t index, const std::vector<double> &F,
                             const std::vector<double> &U,
                             const std::vector<double> &V);

  const Date m_epoch;
  Date m_refTime;
  Date m_curTime;
  double m_resolution;
  std::vector<std::string> m_constituentNames;
  std::vector<size_t> m_constituentIndex;
  std::vector<double> m_latgrid;
  std::vector<std::vector<Tide>> m_tides;
  std::tuple<double, double, double> getStaticParameters(size_t index);

  static constexpr std::array<double, 4> sc{270.434164, 13.1763965268,
                                            -0.0000850, 0.000000039};
  static constexpr std::array<double, 4> hc{279.696678, 0.9856473354,
                                            0.00002267, 0.000000000};
  static constexpr std::array<double, 4> pc{334.329556, 0.1114040803,
                                            -0.0007739, -0.00000026};
  static constexpr std::array<double, 4> npc{-259.183275, 0.0529539222,
                                             -0.0001557, -0.000000050};
  static constexpr std::array<double, 4> ppc{281.220844, 0.0000470684,
                                             0.0000339, 0.000000070};
};

#endif  // TIDEFAC_H
