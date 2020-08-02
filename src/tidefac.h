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
#include <vector>

#include "date.h"
#include "tidefac_global.h"

class TideFac {
 public:
  TIDEFAC_EXPORT TideFac();
  TIDEFAC_EXPORT TideFac(const Date &refTime);
  TIDEFAC_EXPORT TideFac(const double latmin, const double latmax,
                         const double resolution);

  int TIDEFAC_EXPORT generateLatitudeGrid(const double latmin,
                                          const double latmax,
                                          const double resolution);

  int TIDEFAC_EXPORT addConstituent(const std::string &harmonic);
  int TIDEFAC_EXPORT addConstituent(const char *harmonic);
  int TIDEFAC_EXPORT addMajor8();

  void TIDEFAC_EXPORT setRefTime(const Date &refTime);
  Date TIDEFAC_EXPORT refTime() const;

  void TIDEFAC_EXPORT calculate(const Date &d, const double latitude);
  void TIDEFAC_EXPORT calculate(const size_t dt, const double latitude);
  void TIDEFAC_EXPORT calculate(const double dt, const double latitude);
  void TIDEFAC_EXPORT calculate(const Date &d1, const Date &d2,
                                const double latitude);
  void TIDEFAC_EXPORT calculate(const Date &d);
  void TIDEFAC_EXPORT calculate(const size_t dt);
  void TIDEFAC_EXPORT calculate(const double dt);

  void TIDEFAC_EXPORT show(const size_t index = 0) const;

  Date TIDEFAC_EXPORT curTime() const;

  std::string TIDEFAC_EXPORT name(size_t index, size_t gridIndex = 0);
  double TIDEFAC_EXPORT amplitude(size_t index, size_t gridIndex = 0);
  double TIDEFAC_EXPORT frequency(size_t index, size_t gridIndex = 0);
  double TIDEFAC_EXPORT earthTideReductionFactor(size_t index,
                                                 size_t gridIndex = 0);
  double TIDEFAC_EXPORT nodeFactor(size_t index, size_t gridIndex = 0);
  double TIDEFAC_EXPORT equilibriumArgument(size_t index, size_t gridIndex = 0);
  double TIDEFAC_EXPORT nodefactorCorrection(size_t index,
                                             size_t gridIndex = 0);
  double TIDEFAC_EXPORT astronomicArgument(size_t index, size_t gridIndex = 0);

  int TIDEFAC_EXPORT getInterpolationFactors(const double &latitude,
                                             unsigned long &gridIndex,
                                             double &weight);

 private:
  struct Tide {
    Tide(const std::string &name, const double amp, const double freq,
         const double etrf, const double nodefactor, const double eqarg,
         const double nodecorrection, const double astroarg)
        : name(name),
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

  static std::string toUpper(const std::string &s);

  static constexpr double pi() { return 3.14159265358979323846264338327950288; }
  static constexpr double twopi() { return pi() * 2.0; }
  static constexpr double pi180() { return pi() / 180.0; }

  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
  computeAstronomicalArguments(const std::array<double, 6> &astro,
                               const double latitude);

  std::array<double, 162> computeUU(const std::array<double, 6> &astro);

  std::array<std::complex<double>, 162> computeMat(
      const std::array<double, 162> &uu, const std::array<double, 162> &rr);

  std::array<double, 162> computeRR(const double latitude);

  void computeOrbitalParameters(const Date &d, std::array<double, 6> *astro,
                                std::array<double, 6> *ader = nullptr);

  void computePrimaryFactors(std::vector<std::complex<double>> &F,
                             std::vector<double> &U, std::vector<double> &V,
                             const std::array<std::complex<double>, 162> &mat,
                             const std::array<double, 6> &astro);

  void computeMinorCorrections(std::vector<std::complex<double>> &F,
                               std::vector<double> &U, std::vector<double> &V);

  std::vector<double> complexToReal(std::vector<std::complex<double>> &complex);

  std::vector<double> generateGrid(const double min, const double max,
                                   const double res);

  void reorientAngularParameters(double &nc, double &eq, double &aa);

  TideFac::Tide generateTide(const size_t index, const std::vector<double> &F,
                             const std::vector<double> &U,
                             const std::vector<double> &V);

  Date m_refTime;
  Date m_curTime;
  double m_resolution;
  std::vector<std::string> m_constituentNames;
  std::vector<size_t> m_constituentIndex;
  std::vector<double> m_latgrid;
  std::vector<std::vector<Tide>> m_tides;
  std::tuple<double, double, double> getStaticParameters(size_t index);
};

#endif  // TIDEFAC_H
