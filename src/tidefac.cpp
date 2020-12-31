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
#include "tidefac.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <tuple>

#include "constituent.h"

/**
 * @brief Default constructor
 *
 * Initializes to epoch time with a single tide calculation location
 *
 */
TideFac::TideFac()
    : m_refTime(0), m_curTime(0), m_resolution(0.0), m_tides(1) {}

/**
 * @brief Constructor that initializes with a reference time
 * @param refTime initializes the reference time
 */
TideFac::TideFac(const Date &refTime)
    : m_refTime(refTime), m_curTime(0), m_resolution(0.0), m_tides(1) {}

/**
 * @brief Constructor to initialize with a tide grid
 * @param latmin minimum latitude for grid
 * @param latmax maximum latitude for grid
 * @param resolution tide grid resolution
 */
TideFac::TideFac(const double latmin, const double latmax,
                 const double resolution)
    : m_refTime(0),
      m_curTime(0),
      m_resolution(resolution),
      m_latgrid(generateGrid(latmin, latmax, resolution)),
      m_tides(m_latgrid.size()) {}

/**
 * @brief Generates the tide grid with specified coordinates
 * @param latmin minimum latitude for grid
 * @param latmax maximum latitude for grid
 * @param resolution tide grid resolution
 * @return error if not 0
 */
int TideFac::generateLatitudeGrid(const double latmin, const double latmax,
                                  const double resolution) {
  assert(resolution > 0.0);
  assert(latmax > latmin);
  assert(latmin >= -90.0);
  assert(latmax <= 90.0);
  this->m_resolution = resolution;
  this->m_latgrid = TideFac::generateGrid(latmin, latmax, resolution);
  this->m_tides.clear();
  this->m_tides.resize(this->m_latgrid.size());
  return 0;
}

/**
 * @brief Generates the interpolation factors for a given location on the tide
 * grid
 * @param latitude query location
 * @param gridIndex returned lower grid index
 * @param weight weight of the lower grid index
 * @return error if not 0
 */
int TideFac::getInterpolationFactors(const double latitude,
                                     unsigned long &gridIndex,
                                     double &weight) const {
  auto p = latitude < 0.0 ? std::upper_bound(this->m_latgrid.begin(),
                                             this->m_latgrid.end(), latitude)
                          : std::lower_bound(this->m_latgrid.begin(),
                                             this->m_latgrid.end(), latitude);
  gridIndex = p - this->m_latgrid.begin() - 1;
  weight = 1.0 - (latitude - this->m_latgrid[gridIndex]) / this->m_resolution;
  return 0;
}

/**
 * @brief Converts the input string to uppercase
 * @param s string to convert
 * @return uppercase string
 */
void TideFac::toUpper(std::string &s) {
  std::locale loc;
  for (auto &c : s) {
    c = std::toupper(c, loc);
  }
}

/**
 * @brief Add a constituent to the calculation
 * @param harmonic character name of the harmonic
 * @return error if not 0
 */
int TideFac::addConstituent(const char *harmonic) {
  return this->addConstituent(std::string(harmonic));
}

/**
 * @overload
 * @brief Add a constituent to the calculation
 * @param harmonic
 * @return
 */
int TideFac::addConstituent(std::string harmonic) {
  TideFac::toUpper(harmonic);
  size_t idx = Constituent::constituentIndex(harmonic.c_str());
  if (idx == Constituent::null_value<size_t>()) {
    return 1;
  } else {
    auto s = std::find(this->m_constituentIndex.begin(),
                       this->m_constituentIndex.end(), idx);
    if (this->m_constituentIndex.empty() ||
        s == this->m_constituentIndex.end()) {
      this->m_constituentIndex.push_back(idx);
      this->m_constituentNames.push_back(harmonic);
    }
    return 0;
  }
}

/**
 * @brief Shortcut to add the major 8 tide constituents
 * @return error if not 0
 */
int TideFac::addMajor8() {
  int ierr = this->addConstituent("M2");
  ierr += this->addConstituent("S2");
  ierr += this->addConstituent("N2");
  ierr += this->addConstituent("K2");
  ierr += this->addConstituent("K1");
  ierr += this->addConstituent("O1");
  ierr += this->addConstituent("Q1");
  ierr += this->addConstituent("P1");
  return ierr;
}

/**
 * @brief Computes the tides for the given time
 * @param dt time (seconds) since reference time for the calculation
 * @param latitude location for the calculation
 */
void TideFac::calculate(const size_t dt, const double latitude) {
  this->calculate(this->m_refTime + dt, latitude);
}

/**
 * @brief Computes the tides for the given time
 * @overload
 * @param dt time (seconds0 since the reference time
 * @param latitude location for the calculation
 *
 * Note that resolution more fine than one second is not supported
 */
void TideFac::calculate(const double dt, const double latitude) {
  this->calculate(static_cast<size_t>(std::floor(dt)), latitude);
}

std::tuple<double, double, double> TideFac::getStaticParameters(
    const size_t index) {
  return {std::abs(Constituent::constituents()
                       .at(this->m_constituentIndex[index])
                       .doodson_amplitude) *
              0.2675,
          Constituent::constituents()
                  .at(this->m_constituentIndex[index])
                  .frequency *
              TideFac::twopi() / 3600.0,
          Constituent::constituents()
              .at(this->m_constituentIndex[index])
              .earth_tide_reduction_factor};
}

/**
 * @brief Calculation function which uses two dates and returns a set of mean
 * values to use
 * @param d1 start date for simulation
 * @param d2 end date for simulation
 * @param latitude location for the calculation
 */
void TideFac::calculate(const Date &d1, const Date &d2, const double latitude) {
  assert(!this->m_constituentIndex.empty());
  if (this->m_constituentIndex.empty()) {
    std::cerr << "[ERROR]: No tide selected. Aborting calculation" << std::endl;
    return;
  }
  std::vector<Date> dateList;
  size_t dt = std::floor((d2.toSeconds() - d1.toSeconds()) /
                         TideFac::meanNumTides<size_t>());
  dateList.resize(TideFac::meanNumTides<size_t>());
  for (size_t i = 0; i < dateList.size(); ++i) {
    dateList[i] = d1 + i * dt;
  }

  std::vector<double> mean_nodeFactor(this->m_constituentIndex.size(), 0.0);
  std::vector<double> mean_nfCorrection(this->m_constituentIndex.size(), 0.0);
  std::vector<double> V1;

  for (size_t i = 0; i < TideFac::meanNumTides<size_t>(); ++i) {
    std::array<double, 6> astro = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    TideFac::computeOrbitalParameters(dateList[i], &astro);

    std::vector<double> F, U, V;
    TideFac::computeAstronomicalArguments(astro, latitude, F, U, V);

    if (i == 0) V1 = V;  // astro arg at start used

    for (size_t j = 0; j < this->m_constituentIndex.size(); ++j) {
      mean_nodeFactor[j] += F[this->m_constituentIndex[j]];
      mean_nfCorrection[j] += U[this->m_constituentIndex[j]];
    }
  }

  this->m_tides[0].clear();
  this->m_tides.reserve(this->m_constituentIndex.size());
  for (size_t i = 0; i < this->m_constituentIndex.size(); ++i) {
    double nf = mean_nodeFactor[i] / TideFac::meanNumTides<double>();
    double nc =
        (mean_nfCorrection[i] / TideFac::meanNumTides<double>()) * 360.0;
    double eq = (mean_nfCorrection[i] / TideFac::meanNumTides<double>() +
                 V1[this->m_constituentIndex[i]]) *
                360.0;
    double aa = V1[this->m_constituentIndex[i]] * 360.0;
    TideFac::reorientAngularParameters(nc, eq, aa);

    double amp, frequency, etrf;
    std::tie(amp, frequency, etrf) = getStaticParameters(i);

    this->m_tides[0].push_back(Tide(this->m_constituentNames[i], amp, frequency,
                                    etrf, nf, eq, nc, aa));
  }

  this->m_curTime = d1 + d2.toSeconds() / 2;
}

/**
 * @overload
 * @brief Computes the values for the given time at all locations on the
 * latitude grid
 * @param d date for the calculation
 */
void TideFac::calculate(const Date &d) {
  assert(!this->m_constituentIndex.empty());
  if (this->m_constituentIndex.empty()) {
    std::cerr << "[ERROR]: No tide selected. Aborting calculation" << std::endl;
    return;
  }

  assert(this->m_tides.size() > 1);
  if (this->m_tides.size() <= 1) {
    std::cerr << "[ERROR] Calculation not configured for tide grid."
              << std::endl;
    return;
  }

  this->m_curTime = d;

  std::array<double, 6> astro = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  TideFac::computeOrbitalParameters(d, &astro);

  for (size_t i = 0; i < this->m_tides.size(); ++i) {
    std::vector<double> F, U, V;
    TideFac::computeAstronomicalArguments(astro, this->m_latgrid[i], F, U, V);

    this->m_tides[i].clear();
    this->m_tides[i].reserve(this->m_constituentIndex.size());
    for (size_t j = 0; j < this->m_constituentIndex.size(); ++j) {
      this->m_tides[i].push_back(this->generateTide(j, F, U, V));
    }
  }
}

/**
 * @overload
 * @brief Computes all points on the tide grid for a specified time since the
 * reference time
 * @param dt time (seconds) since the reference
 */
void TideFac::calculate(const size_t dt) {
  return this->calculate(this->m_refTime + dt);
}

/**
 * @overload
 * @brief Computes all points on the tide grid for a specified time since the
 * reference time
 * @param dt time (seconds) since the reference
 */
void TideFac::calculate(const double dt) {
  this->calculate(static_cast<size_t>(std::floor(dt)));
}

/**
 * @overload
 * @brief Computes the tide parameters for a given date
 * @param d date for the calculation
 * @param latitude location for the calculation
 */
void TideFac::calculate(const Date &d, const double latitude) {
  assert(!this->m_constituentIndex.empty());
  if (this->m_constituentIndex.empty()) {
    std::cerr << "[INFO]: No tide selected. Aborting calculation" << std::endl;
    return;
  }

  std::array<double, 6> astro = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  TideFac::computeOrbitalParameters(d, &astro);

  std::vector<double> F, U, V;
  TideFac::computeAstronomicalArguments(astro, latitude, F, U, V);

  this->m_curTime = d;

  this->m_tides[0].clear();
  this->m_tides.reserve(this->m_constituentIndex.size());
  for (size_t i = 0; i < this->m_constituentIndex.size(); ++i) {
    this->m_tides[0].push_back(this->generateTide(i, F, U, V));
  }
}

/**
 * @brief Generates the tide calculation for the given F,U,V
 * @param index tide index
 * @param F tide factor vector
 * @param U tide factor vector
 * @param V tide factor vector
 * @return Tide object for the specified index
 */
TideFac::Tide TideFac::generateTide(const size_t index,
                                    const std::vector<double> &F,
                                    const std::vector<double> &U,
                                    const std::vector<double> &V) {
  double nodeFactor = F[this->m_constituentIndex[index]];
  double nc = U[this->m_constituentIndex[index]] * 360.0;
  double eq = (U[this->m_constituentIndex[index]] +
               V[this->m_constituentIndex[index]]) *
              360.0;
  double aa = V[this->m_constituentIndex[index]] * 360.0;

  TideFac::reorientAngularParameters(nc, eq, aa);

  double amp, frequency, etrf;
  std::tie(amp, frequency, etrf) = this->getStaticParameters(index);

  return Tide(this->m_constituentNames[index], amp, frequency, etrf, nodeFactor,
              eq, nc, aa);
}

/**
 * @brief Compute the astronomic arguments for the given orbital arguments
 * @param astro orbital parameters
 * @param latitude location for the calculation
 * @return three vectors of tide factors in a tuple
 */
void TideFac::computeAstronomicalArguments(const std::array<double, 6> &astro,
                                           const double latitude,
                                           std::vector<double> &F,
                                           std::vector<double> &U,
                                           std::vector<double> &V) {
  std::array<double, 162> rr = TideFac::computeRR(latitude);
  std::array<double, 162> uu = TideFac::computeUU(astro);
  std::array<std::complex<double>, 162> mat = TideFac::computeMat(uu, rr);

  std::vector<std::complex<double>> Fcomplex(
      Constituent::constituents().size());
  F.resize(Constituent::constituents().size());
  U.resize(Constituent::constituents().size());
  V.resize(Constituent::constituents().size());

  TideFac::computePrimaryFactors(Fcomplex, U, V, mat, astro);
  TideFac::computeMinorCorrections(Fcomplex, U, V);
  F = TideFac::complexToReal(Fcomplex);
}

/**
 * @brief Computes the orbital parameters for the given date
 * @param d time for the calculation
 * @param astro astronomic arguments
 * @param ader astronomic arguments
 */
void TideFac::computeOrbitalParameters(const Date &d,
                                       std::array<double, 6> *astro,
                                       std::array<double, 6> *ader) {
  const Date epoch(1899, 12, 31, 12, 0, 0);
  constexpr std::array<double, 4> sc{270.434164, 13.1763965268, -0.0000850,
                                     0.000000039};
  constexpr std::array<double, 4> hc{279.696678, 0.9856473354, 0.00002267,
                                     0.000000000};
  constexpr std::array<double, 4> pc{334.329556, 0.1114040803, -0.0007739,
                                     -0.00000026};
  constexpr std::array<double, 4> npc{-259.183275, 0.0529539222, -0.0001557,
                                      -0.000000050};
  constexpr std::array<double, 4> ppc{281.220844, 0.0000470684, 0.0000339,
                                      0.000000070};

  const Date jd(d - epoch.toSeconds());
  const double dd = static_cast<double>(jd.toSeconds()) / 86400.0;
  const double nd = dd / 10000.0;
  const std::array<double, 4> args{1.0, dd, nd * nd, nd * nd * nd};

  astro->at(1) = std::fmod(
      (sc[0] * args[0] + sc[1] * args[1] + sc[2] * args[2] + sc[3] * args[3]) /
          360.0,
      1.0);
  astro->at(2) = std::fmod(
      (hc[0] * args[0] + hc[1] * args[1] + hc[2] * args[2] + hc[3] * args[3]) /
          360.0,
      1.0);
  astro->at(3) = std::fmod(
      (pc[0] * args[0] + pc[1] * args[1] + pc[2] * args[2] + pc[3] * args[3]) /
          360.0,
      1.0);
  astro->at(4) = std::fmod((npc[0] * args[0] + npc[1] * args[1] +
                            npc[2] * args[2] + npc[3] * args[3]) /
                               360.0,
                           1.0);
  astro->at(5) = std::fmod((ppc[0] * args[0] + ppc[1] * args[1] +
                            ppc[2] * args[2] + ppc[3] * args[3]) /
                               360.0,
                           1.0);
  astro->at(0) = std::fmod(static_cast<double>(d.toSeconds()) / 86400.0, 1.0) +
                 astro->at(2) - astro->at(1);

  if (ader == nullptr) return;

  const std::array<double, 4> args2 = {0.0, 1.0, 2.0e-4 * nd, 3.0e-4 * nd * nd};
  ader->at(1) = std::fmod((sc[0] * args[0] + sc[1] * args2[1] +
                           sc[2] * args2[2] + sc[3] * args2[3]) /
                              360.0,
                          1.0);
  ader->at(2) = std::fmod((hc[0] * args2[0] + hc[1] * args2[1] +
                           hc[2] * args2[2] + hc[3] * args2[3]) /
                              360.0,
                          1.0);
  ader->at(3) = std::fmod((pc[0] * args2[0] + pc[1] * args2[1] +
                           pc[2] * args2[2] + pc[3] * args2[3]) /
                              360.0,
                          1.0);
  ader->at(4) = std::fmod((npc[0] * args2[0] + npc[1] * args2[1] +
                           npc[2] * args2[2] + npc[3] * args2[3]) /
                              360.0,
                          1.0);
  ader->at(5) = std::fmod((ppc[0] * args2[0] + ppc[1] * args2[1] +
                           ppc[2] * args2[2] + ppc[3] * args2[3]) /
                              360.0,
                          1.0);
  ader->at(0) = 1.0 + ader->at(2) - ader->at(1);
}

/**
 * @brief Quickly prints the specified tide on screen
 * @param index index for the tide object to show
 */
void TideFac::show(const size_t index) const {
  std::cout << "Tidal Factors for simulation time " << this->m_curTime << "\n";
  for (auto &t : this->m_tides) {
    std::cout << t[index].name << " " << t[index].amp << " " << t[index].freq
              << " " << t[index].etrf << " " << t[index].nodefactor << " "
              << t[index].eqarg << "\n";
  }
}

/**
 * @brief Computes the matrix sum used in the calculation where constituent
 * indicies align
 * @param mat input matrix
 * @param idx index of the tide
 * @return sum
 */
std::complex<double> matsum(const std::array<std::complex<double>, 162> &mat,
                            const int idx) {
  std::complex<double> s = {0.0, 0.0};
  for (size_t i = 0; i < mat.size(); ++i) {
    if (Constituent::iconst().at(i) == idx) {
      s += mat[i];
    }
  }
  return s;
}

/**
 * @brief Computes the U components
 * @param astro astronomic arguments
 * @return array of U for each tide
 */
std::array<double, 162> TideFac::computeUU(const std::array<double, 6> &astro) {
  std::array<double, 162> uu{};
  uu.fill(0.0);

  for (size_t i = 0; i < uu.size(); ++i) {
    uu[i] = std::fmod((Constituent::deldood().at(i)[0] * astro[3] +
                       Constituent::deldood().at(i)[1] * astro[4] +
                       Constituent::deldood().at(i)[2] * astro[5]) +
                          Constituent::phcorr().at(i),
                      1.0);
  }
  return uu;
}

/**
 * @brief Computes the matrix rr*e^(i*pi*uu)
 * @param uu U array
 * @param rr R array
 * @return matrix containing values to be summed later
 */
std::array<std::complex<double>, 162> TideFac::computeMat(
    const std::array<double, 162> &uu, const std::array<double, 162> &rr) {
  constexpr std::complex<double> c_i(0, 1);
  std::array<std::complex<double>, 162> mat;
  for (size_t i = 0; i < rr.size(); ++i) {
    mat[i] = rr[i] * std::exp(c_i * TideFac::twopi() * uu[i]);
  }
  return mat;
}

/**
 * @brief Compute the RR array
 * @param latitude location for calculation
 * @return RR array
 */
std::array<double, 162> TideFac::computeRR(const double latitude) {
  auto latsign = [](double lat) {
    return std::abs(lat) < 0.5 ? std::copysign(0.5, lat) : lat;
  };

  const double slat = std::sin(TideFac::deg2rad() * latsign(latitude));
  const double slatfac1 = 0.36309 * (1.0 - 5.0 * slat * slat) / slat;
  const double slatfac2 = slat * 2.59808;

  std::array<double, 162> rr = Constituent::amprat();
  for (size_t i = 0; i < Constituent::amprat().size(); ++i) {
    if (Constituent::ilatfac().at(i) == 1) {
      rr[i] = Constituent::amprat().at(i) * slatfac1;
    } else if (Constituent::ilatfac().at(i) == 2) {
      rr[i] = rr[i] * slatfac2;
    }
  }
  return rr;
}

/**
 * @brief Computes the minor corrections for F,U,V vectors in place
 */
void TideFac::computeMinorCorrections(std::vector<std::complex<double>> &F,
                                      std::vector<double> &U,
                                      std::vector<double> &V) {
  for (size_t i = 0; i < Constituent::constituents().size(); ++i) {
    if (Constituent::constituents().at(i).i_shallow !=
        Constituent::null_value<int>()) {
      size_t j = Constituent::constituents().at(i).n_shallow;
      size_t k = Constituent::constituents().at(i).i_shallow;
      std::vector<int> ik(j), iloc(j);
      std::vector<double> exp1(j), exp2(j);
      for (size_t p = 0; p < ik.size(); ++p) {
        ik[p] = k + p;
        iloc[p] = Constituent::shallow_iname().at(ik[p] - 1);
        exp1[p] = Constituent::shallow_coef().at(ik[p] - 1);
        exp2[p] = std::abs(exp1[p]);
      }

      double product_1 = 1.0;
      double product_2 = 0.0;
      double product_3 = 0.0;
      for (size_t p = 0; p < ik.size(); ++p) {
        product_1 = product_1 * std::pow(F[iloc[p] - 1].real(), exp2[p]);
        product_2 = product_2 + U[iloc[p] - 1] * exp1[p];
        product_3 = product_3 + V[iloc[p] - 1] * exp1[p];
      }
      F[i] = product_1;
      U[i] = product_2;
      V[i] = product_3;
    }
  }
}

/**
 * @brief Computes the primary values for F, U, V vectors
 */
void TideFac::computePrimaryFactors(
    std::vector<std::complex<double>> &F, std::vector<double> &U,
    std::vector<double> &V, const std::array<std::complex<double>, 162> &mat,
    const std::array<double, 6> &astro) {
  auto ind = Constituent::iconst_unique();
  for (int i : ind) {
    F[i - 1] = 1.0 + matsum(mat, i);
  }

  for (size_t i = 0; i < Constituent::constituents().size(); ++i) {
    U[i] = std::log(F[i]).imag() / TideFac::twopi();
    F[i] = std::abs(F[i]);
    if (Constituent::constituents().at(i).doodson[0] !=
        Constituent::null_value<double>()) {
      V[i] = std::fmod(
          Constituent::constituents().at(i).doodson[0] * astro[0] +
              Constituent::constituents().at(i).doodson[1] * astro[1] +
              Constituent::constituents().at(i).doodson[2] * astro[2] +
              Constituent::constituents().at(i).doodson[3] * astro[3] +
              Constituent::constituents().at(i).doodson[4] * astro[4] +
              Constituent::constituents().at(i).doodson[5] * astro[5] +
              Constituent::constituents().at(i).semi,
          1.0);
    } else {
      V[i] = Constituent::null_value<double>();
    }
  }
}

/**
 * @brief Converts a complex double to a double
 * @param complex vector containing complex values
 * @return vector containing double values
 */
std::vector<double> TideFac::complexToReal(
    std::vector<std::complex<double>> &complex) {
  std::vector<double> dbl(complex.size());
  std::transform(complex.begin(), complex.end(), dbl.begin(),
                 [](std::complex<double> &d) { return d.real(); });
  return dbl;
}

/**
 * @brief Generates the latitude grid based on min, max, res
 * @param min minimum value in the grid
 * @param max maximum value in the grid
 * @param res resolution of grid
 * @return vector containing grid points
 */
std::vector<double> TideFac::generateGrid(const double min, const double max,
                                          const double res) {
  size_t n = std::floor((max - min) / res) + 2;
  std::vector<double> l;
  l.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    l.push_back(min + (i - 1) * res);
  }
  return l;
}

/**
 * @brief Converts an angle to 360 degree orientation
 * @param v value to be converted in place
 */
inline void zeroTo360(double &v) {
  if (v < 0.0) v += 360.0;
}

/**
 * @brief Converts the angles in the calculation to 0-360 degrees in place
 * @param nc node correction
 * @param eq equilibrium argument
 * @param aa astronomic argument
 */
void TideFac::reorientAngularParameters(double &nc, double &eq, double &aa) {
  zeroTo360(nc);
  zeroTo360(eq);
  zeroTo360(aa);
}

/**
 * @brief Returns the tide name at the specified position in the grid
 * @param index tide index
 * @param gridIndex grid index
 * @return name of tide
 */
std::string TideFac::name(size_t index, size_t gridIndex) const {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].name;
}

/**
 * @brief Amplitude of tide at the specified position
 * @param index tide index
 * @param gridIndex grid index
 * @return amplitude of tide
 */
double TideFac::amplitude(size_t index, size_t gridIndex) const {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].amp;
}

/**
 * @brief Returns the frequency of the tide at the specified position in the
 * grid
 * @param index tide index
 * @param gridIndex grid index
 * @return frequency of tide
 */
double TideFac::frequency(size_t index, size_t gridIndex) const {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].freq;
}

/**
 * @brief Returns the earth tide reduction factor at the specified position in
 * the grid
 * @param index tide index
 * @param gridIndex grid index
 * @return earth tide reduction factor for the tide
 */
double TideFac::earthTideReductionFactor(size_t index, size_t gridIndex) const {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].etrf;
}

/**
 * @brief Returns the node factor of the tide at the specified position in the
 * grid
 * @param index tide index
 * @param gridIndex grid index
 * @return node factor of tide
 */
double TideFac::nodeFactor(size_t index, size_t gridIndex) const {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].nodefactor;
}

/**
 * @brief Returns the equilibrium argument of the tide at the specified position
 * in the grid
 * @param index tide index
 * @param gridIndex grid index
 * @return equilibrium argument of tide
 */
double TideFac::equilibriumArgument(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].eqarg;
}

/**
 * @brief Returns the node factor correction of the tide at the specified
 * position in the grid
 * @param index tide index
 * @param gridIndex grid index
 * @return node factor correction for the tide
 */
double TideFac::nodefactorCorrection(size_t index, size_t gridIndex) const {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].nodecorrection;
}

/**
 * @brief Returns the astronomic argument of the tide at the specified position
 * in the grid
 * @param index tide index
 * @param gridIndex grid index
 * @return astronomic argument of tide
 */
double TideFac::astronomicArgument(size_t index, size_t gridIndex) const {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].astroarg;
}

/**
 * @brief Reference time for the current calculation
 * @return reference time
 */
Date TideFac::refTime() const { return this->m_refTime; }

/**
 * @brief Sets the reference time for the calculation
 * @param refTime date object with reference time
 */
void TideFac::setRefTime(const Date &refTime) { this->m_refTime = refTime; }

/**
 * @brief Returns the current time for the calculation
 * @return date object with current time
 */
Date TideFac::curTime() const { return this->m_curTime; }
