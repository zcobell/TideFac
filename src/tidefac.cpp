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
#include <locale>
#include <numeric>

#include "constituent.h"

TideFac::TideFac()
    : m_refTime(0), m_curTime(0), m_resolution(0.0), m_tides(1) {}

TideFac::TideFac(const Date &refTime)
    : m_refTime(refTime), m_curTime(0), m_resolution(0.0), m_tides(1) {}

TideFac::TideFac(const double latmin, const double latmax,
                 const double resolution)
    : m_refTime(0),
      m_curTime(0),
      m_resolution(resolution),
      m_latgrid(generateGrid(latmin, latmax, resolution)),
      m_tides(m_latgrid.size()) {}

int TideFac::generateLatitudeGrid(const double latmin, const double latmax,
                                  const double resolution) {
  assert(resolution > 0.0);
  assert(latmax > latmin);
  assert(latmin >= -90.0);
  assert(latmax <= 90.0);
  this->m_resolution = resolution;
  this->m_latgrid = this->generateGrid(latmin, latmax, resolution);
  this->m_tides.clear();
  this->m_tides.resize(this->m_latgrid.size());
  return 0;
}

int TideFac::addConstituent(const char *harmonic) {
  return this->addConstituent(std::string(harmonic));
}

int TideFac::getInterpolationFactors(const double &latitude,
                                     unsigned long &gridIndex, double &weight) {
  if (latitude < 0.0) {
    auto p = std::upper_bound(this->m_latgrid.begin(), this->m_latgrid.end(),
                              latitude);
    if (p == this->m_latgrid.begin() || p == this->m_latgrid.end()) return 1;
    gridIndex = p - this->m_latgrid.begin() - 1;
  } else {
    auto p = std::lower_bound(this->m_latgrid.begin(), this->m_latgrid.end(),
                              latitude);
    if (p == this->m_latgrid.begin() || p == this->m_latgrid.end()) return 1;
    gridIndex = p - this->m_latgrid.begin() - 1;
  }
  weight = 1.0 - (latitude - this->m_latgrid[gridIndex]) / this->m_resolution;
  return 0;
}

std::string TideFac::toUpper(const std::string &s) {
  std::locale loc;
  std::string sout;
  sout.reserve(s.size());
  for (auto &c : s) {
    sout.push_back(std::toupper(c, loc));
  }
  return sout;
}

int TideFac::addConstituent(const std::string &harmonic) {
  const std::string harm = TideFac::toUpper(harmonic);
  size_t idx = Constituent::constituentIndex(harm);
  if (idx == Constituent::nullvalue<size_t>()) {
    return 1;
  } else {
    auto s = std::find(this->m_constituentIndex.begin(),
                       this->m_constituentIndex.end(), idx);
    if (this->m_constituentIndex.size() == 0 ||
        s == this->m_constituentIndex.end()) {
      this->m_constituentIndex.push_back(idx);
      this->m_constituentNames.push_back(harm);
    }
    return 0;
  }
}

int TideFac::addMajor8() {
  this->addConstituent("M2");
  this->addConstituent("S2");
  this->addConstituent("N2");
  this->addConstituent("K2");
  this->addConstituent("K1");
  this->addConstituent("O1");
  this->addConstituent("Q1");
  this->addConstituent("P1");
  return 0;
}

void TideFac::calculate(const size_t dt, const double latitude) {
  this->calculate(this->m_refTime + dt, latitude);
}

void TideFac::calculate(const double dt, const double latitude) {
  this->calculate(static_cast<size_t>(std::floor(dt)), latitude);
}

std::tuple<double, double, double> TideFac::getStaticParameters(
    const size_t index) {
  return {std::abs(Constituent::constituents()
                       ->at(this->m_constituentIndex[index])
                       ->doodsonamp) *
              0.2675,
          Constituent::constituents()
                  ->at(this->m_constituentIndex[index])
                  ->frequency *
              TideFac::twopi() / 3600.0,
          Constituent::constituents()
              ->at(this->m_constituentIndex[index])
              ->earthreduc};
}

void zeroTo360(double &v) {
  if (v < 0.0) v += 360.0;
}

void TideFac::reorientAngularParameters(double &nc, double &eq, double &aa) {
  zeroTo360(nc);
  zeroTo360(eq);
  zeroTo360(aa);
}

void TideFac::calculate(const Date &d1, const Date &d2, const double latitude) {
  if (this->m_constituentIndex.size() == 0) {
    std::cerr << "[ERROR]: No tide selected. Aborting calculation" << std::endl;
    return;
  }
  std::vector<Date> dateList;
  size_t dt = std::floor((d2.toSeconds() - d1.toSeconds()) / 1000.0);
  dateList.resize(1000);
  for (size_t i = 0; i < dateList.size(); ++i) {
    dateList[i] = d1 + i * dt;
  }

  std::vector<double> mean_nodeFactor(this->m_constituentIndex.size(), 0.0);
  std::vector<double> mean_nfCorrection(this->m_constituentIndex.size(), 0.0);
  std::vector<double> V1;

  for (size_t i = 0; i < 1000; ++i) {
    std::array<double, 6> astro;
    this->computeOrbitalParameters(dateList[i], &astro);

    std::vector<double> F, U, V;
    std::tie(F, U, V) = this->computeAstronomicalArguments(astro, latitude);

    if (i == 0) V1 = V;  // astro arg at start used

    for (size_t j = 0; j < this->m_constituentIndex.size(); ++j) {
      mean_nodeFactor[j] += F[this->m_constituentIndex[j]];
      mean_nfCorrection[j] += U[this->m_constituentIndex[j]];
    }
  }

  this->m_tides[0].clear();
  this->m_tides.reserve(this->m_constituentIndex.size());
  for (size_t i = 0; i < this->m_constituentIndex.size(); ++i) {
    double nf = mean_nodeFactor[i] / 1000.0;
    double nc = (mean_nfCorrection[i] / 1000.0) * 360.0;
    double eq =
        (mean_nfCorrection[i] / 1000.0 + V1[this->m_constituentIndex[i]]) *
        360.0;
    double aa = V1[this->m_constituentIndex[i]] * 360.0;
    this->reorientAngularParameters(nc, eq, aa);

    double amp, frequency, etrf;
    std::tie(amp, frequency, etrf) = getStaticParameters(i);

    this->m_tides[0].push_back(Tide(this->m_constituentNames[i], amp, frequency,
                                    etrf, nf, eq, nc, aa));
  }

  this->m_curTime = d1 + d2.toSeconds() / 2;
}

void TideFac::calculate(const Date &d) {
  if (this->m_constituentIndex.size() == 0) {
    std::cerr << "[ERROR]: No tide selected. Aborting calculation" << std::endl;
    return;
  }
  if (this->m_tides.size() <= 1) {
    std::cerr << "[ERROR] Calculation not configured for tide grid."
              << std::endl;
    return;
  }

  this->m_curTime = d;

  std::array<double, 6> astro;
  this->computeOrbitalParameters(d, &astro);

  for (size_t i = 0; i < this->m_tides.size(); ++i) {
    std::vector<double> F, U, V;
    std::tie(F, U, V) =
        this->computeAstronomicalArguments(astro, this->m_latgrid[i]);

    this->m_tides[i].clear();
    this->m_tides[i].reserve(this->m_constituentIndex.size());
    for (size_t j = 0; j < this->m_constituentIndex.size(); ++j) {
      this->m_tides[i].push_back(this->generateTide(j, F, U, V));
    }
  }
}

void TideFac::calculate(const size_t dt) {
  return this->calculate(this->m_refTime + dt);
}

void TideFac::calculate(const double dt) {
  this->calculate(static_cast<size_t>(std::floor(dt)));
}

void TideFac::calculate(const Date &d, const double latitude) {
  if (this->m_constituentIndex.size() == 0) {
    std::cerr << "[INFO]: No tide selected. Aborting calculation" << std::endl;
    return;
  }

  std::array<double, 6> astro;
  this->computeOrbitalParameters(d, &astro);

  std::vector<double> F, U, V;
  std::tie(F, U, V) = this->computeAstronomicalArguments(astro, latitude);

  this->m_curTime = d;

  this->m_tides[0].clear();
  this->m_tides.reserve(this->m_constituentIndex.size());
  for (size_t i = 0; i < this->m_constituentIndex.size(); ++i) {
    this->m_tides[0].push_back(this->generateTide(i, F, U, V));
  }
}

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

  this->reorientAngularParameters(nc, eq, aa);

  double amp, frequency, etrf;
  std::tie(amp, frequency, etrf) = this->getStaticParameters(index);

  return Tide(this->m_constituentNames[index], amp, frequency, etrf, nodeFactor,
              eq, nc, aa);
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
TideFac::computeAstronomicalArguments(const std::array<double, 6> &astro,
                                      const double latitude) {
  std::array<double, 162> rr = this->computeRR(latitude);
  std::array<double, 162> uu = this->computeUU(astro);
  std::array<std::complex<double>, 162> mat = this->computeMat(uu, rr);

  std::vector<std::complex<double>> Fcomplex(
      Constituent::constituents()->size());
  std::vector<double> U(Constituent::constituents()->size());
  std::vector<double> V(Constituent::constituents()->size());
  this->computePrimaryFactors(Fcomplex, U, V, mat, astro);
  this->computeMinorCorrections(Fcomplex, U, V);
  std::vector<double> F = TideFac::complexToReal(Fcomplex);

  return {F, U, V};
}

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

  Date jd(d);
  jd.addSeconds(-epoch.toSeconds());
  const double dd = jd.toSeconds() / 86400.0;
  const double nd = dd / 10000.0;
  std::array<double, 4> args{1.0, dd, nd * nd, nd * nd * nd};

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
  astro->at(0) =
      std::fmod(d.toSeconds() / 86400.0, 1.0) + astro->at(2) - astro->at(1);

  if (ader == nullptr) return;

  args = {0.0, 1.0, 2.0e-4 * nd, 3.0e-4 * nd * nd};
  ader->at(1) = std::fmod(
      (sc[0] * args[0] + sc[1] * args[1] + sc[2] * args[2] + sc[3] * args[3]) /
          360.0,
      1.0);
  ader->at(2) = std::fmod(
      (hc[0] * args[0] + hc[1] * args[1] + hc[2] * args[2] + hc[3] * args[3]) /
          360.0,
      1.0);
  ader->at(3) = std::fmod(
      (pc[0] * args[0] + pc[1] * args[1] + pc[2] * args[2] + pc[3] * args[3]) /
          360.0,
      1.0);
  ader->at(4) = std::fmod((npc[0] * args[0] + npc[1] * args[1] +
                           npc[2] * args[2] + npc[3] * args[3]) /
                              360.0,
                          1.0);
  ader->at(5) = std::fmod((ppc[0] * args[0] + ppc[1] * args[1] +
                           ppc[2] * args[2] + ppc[3] * args[3]) /
                              360.0,
                          1.0);
  ader->at(0) = 1.0 + ader->at(2) - ader->at(1);

  return;
}

void TideFac::show(const size_t index) const {
  std::cout << "Tidal Factors for simulation time " << this->m_curTime << "\n";
  for (auto &t : this->m_tides) {
    std::cout << t[index].name << " " << t[index].amp << " " << t[index].freq
              << " " << t[index].etrf << " " << t[index].nodefactor << " "
              << t[index].eqarg << "\n";
  }
}

std::complex<double> matsum(const std::array<std::complex<double>, 162> &mat,
                            const int idx) {
  std::complex<double> s = {0.0, 0.0};
  for (size_t i = 0; i < mat.size(); ++i) {
    if (Constituent::iconst()->at(i) == idx) {
      s += mat[i];
    }
  }
  return s;
}

std::array<double, 162> TideFac::computeUU(const std::array<double, 6> &astro) {
  std::array<double, 162> uu;
  for (size_t i = 0; i < uu.size(); ++i) {
    uu[i] = std::fmod((Constituent::deldood()->at(i)[0] * astro[3] +
                       Constituent::deldood()->at(i)[1] * astro[4] +
                       Constituent::deldood()->at(i)[2] * astro[5]) +
                          Constituent::phcorr()->at(i),
                      1.0);
  }
  return uu;
}

std::array<std::complex<double>, 162> TideFac::computeMat(
    const std::array<double, 162> &uu, const std::array<double, 162> &rr) {
  constexpr std::complex<double> c_i(0, 1);
  std::array<std::complex<double>, 162> mat;
  for (size_t i = 0; i < rr.size(); ++i) {
    mat[i] = rr[i] * std::exp(c_i * TideFac::twopi() * uu[i]);
  }
  return mat;
}

std::array<double, 162> TideFac::computeRR(const double latitude) {
  double lat = latitude;
  if (std::abs(latitude) < 0.5) {
    lat = std::copysign(0.5, latitude);
  }
  double slat = std::sin(TideFac::pi180() * lat);
  std::array<double, 162> rr = *(Constituent::amprat());
  for (size_t i = 0; i < Constituent::amprat()->size(); ++i) {
    if (Constituent::ilatfac()->at(i) == 1) {
      rr[i] = (Constituent::amprat()->at(i) * 0.36309 *
               (1.0 - 5.0 * slat * slat) / slat);
    } else if (Constituent::ilatfac()->at(i) == 2) {
      rr[i] = rr[i] * 2.59808 * slat;
    }
  }
  return rr;
}

void TideFac::computeMinorCorrections(std::vector<std::complex<double>> &F,
                                      std::vector<double> &U,
                                      std::vector<double> &V) {
  for (size_t i = 0; i < Constituent::constituents()->size(); ++i) {
    if (Constituent::constituents()->at(i)->ishallow !=
        Constituent::nullvalue<int>()) {
      int j = Constituent::constituents()->at(i)->nshallow;
      int k = Constituent::constituents()->at(i)->ishallow;
      std::vector<int> ik(j), iloc(j);
      std::vector<double> exp1(j), exp2(j);
      for (size_t p = 0; p < ik.size(); ++p) {
        ik[p] = k + p;
        iloc[p] = Constituent::shallow_iname()->at(ik[p] - 1);
        exp1[p] = Constituent::shallow_coef()->at(ik[p] - 1);
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

void TideFac::computePrimaryFactors(
    std::vector<std::complex<double>> &F, std::vector<double> &U,
    std::vector<double> &V, const std::array<std::complex<double>, 162> &mat,
    const std::array<double, 6> &astro) {
  auto ind = Constituent::iconst_unique();
  for (size_t i = 0; i < ind->size(); ++i) {
    F[ind->at(i) - 1] = 1.0 + matsum(mat, ind->at(i));
  }

  for (size_t i = 0; i < Constituent::constituents()->size(); ++i) {
    U[i] = std::log(F[i]).imag() / TideFac::twopi();
    F[i] = std::abs(F[i]);
    if (Constituent::constituents()->at(i)->doodson[0] !=
        Constituent::nullvalue<double>()) {
      V[i] = std::fmod(
          Constituent::constituents()->at(i)->doodson[0] * astro[0] +
              Constituent::constituents()->at(i)->doodson[1] * astro[1] +
              Constituent::constituents()->at(i)->doodson[2] * astro[2] +
              Constituent::constituents()->at(i)->doodson[3] * astro[3] +
              Constituent::constituents()->at(i)->doodson[4] * astro[4] +
              Constituent::constituents()->at(i)->doodson[5] * astro[5] +
              Constituent::constituents()->at(i)->semi,
          1.0);
    } else {
      V[i] = Constituent::nullvalue<double>();
    }
  }
}

std::vector<double> TideFac::complexToReal(
    std::vector<std::complex<double>> &complex) {
  std::vector<double> dbl(complex.size());
  std::transform(complex.begin(), complex.end(), dbl.begin(),
                 [](std::complex<double> &d) { return d.real(); });
  return dbl;
}

Date TideFac::curTime() const { return this->m_curTime; }

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

std::string TideFac::name(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].name;
}

double TideFac::amplitude(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].amp;
}

double TideFac::frequency(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].freq;
}

double TideFac::earthTideReductionFactor(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].etrf;
}

double TideFac::nodeFactor(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].nodefactor;
}

double TideFac::equilibriumArgument(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].eqarg;
}

double TideFac::nodefactorCorrection(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].nodecorrection;
}

double TideFac::astronomicArgument(size_t index, size_t gridIndex) {
  assert(gridIndex < this->m_tides.size());
  assert(index < this->m_tides[gridIndex].size());
  return this->m_tides[gridIndex][index].astroarg;
}

Date TideFac::refTime() const { return this->m_refTime; }

void TideFac::setRefTime(const Date &refTime) { this->m_refTime = refTime; }
