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
#include <numeric>

#include "constituent.h"

TideFac::TideFac() : m_refTime(0) {}

TideFac::TideFac(const Date &refTime) : m_refTime(refTime) {}

int TideFac::addConstituent(const char *harmonic) {
  return this->addConstituent(std::string(harmonic));
}

int TideFac::addConstituent(const std::string &harmonic) {
  size_t idx = Constituent::constituentIndex(harmonic);
  if (idx == Constituent::nullvalue<size_t>()) {
    return 1;
  } else {
    auto s = std::find(this->m_constituentIndex.begin(),
                       this->m_constituentIndex.end(), idx);
    if (this->m_constituentIndex.size() == 0 ||
        s == this->m_constituentIndex.end()) {
      this->m_constituentIndex.push_back(idx);
      this->m_constituentNames.push_back(harmonic);
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

void TideFac::calculate(const Date &d, const double latitude) {
  if (this->m_constituentIndex.size() == 0) {
    std::cerr << "[INFO]: No tide selected. Aborting calculation" << std::endl;
    return;
  }

  std::vector<double> F, U, V;
  std::tie(F, U, V) = this->computeAstronomicalArguments(d, latitude);

  this->m_curTime = d;

  this->m_tides.clear();
  this->m_tides.reserve(this->m_constituentIndex.size());
  for (size_t i = 0; i < this->m_constituentIndex.size(); ++i) {
    std::string name = this->m_constituentNames[i];
    double nodeFactor = F[this->m_constituentIndex[i]];
    double nfCorrection =
        (U[this->m_constituentIndex[i]] + V[this->m_constituentIndex[i]]) *
        360.0;
    nfCorrection = nfCorrection < 0.0 ? nfCorrection + 360.0 : nfCorrection;
    double frequency = Constituent::constituents()
                           ->at(this->m_constituentIndex[i])
                           ->frequency *
                       TideFac::twopi() / 3600.0;
    double etrf = Constituent::constituents()
                      ->at(this->m_constituentIndex[i])
                      ->earthreduc;
    double amp = std::abs(Constituent::constituents()
                              ->at(this->m_constituentIndex[i])
                              ->doodsonamp) *
                 0.2675;
    this->m_tides.push_back(
        Tide{name, amp, frequency, etrf, nodeFactor, nfCorrection});
  }
}

void TideFac::show() const {
  std::cout << "Tidal Factors for simulation time " << this->m_curTime << "\n";
  for (auto &t : this->m_tides) {
    std::cout << t.name << " " << t.amp << " " << t.freq << " " << t.etrf << " "
              << t.nodefactor << " " << t.eqarg << "\n";
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

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
TideFac::computeAstronomicalArguments(const Date &d, const double latitude) {
  std::array<double, 6> astro;
  this->computeOrbitalParameters(d, &astro);

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

  std::array<double, 162> uu;
  for (size_t i = 0; i < uu.size(); ++i) {
    uu[i] = std::fmod((Constituent::deldood()->at(i)[0] * astro[3] +
                       Constituent::deldood()->at(i)[1] * astro[4] +
                       Constituent::deldood()->at(i)[2] * astro[5]) +
                          Constituent::phcorr()->at(i),
                      1.0);
  }

  constexpr std::complex<double> c_i(0, 1);
  std::array<std::complex<double>, 162> mat;
  for (size_t i = 0; i < rr.size(); ++i) {
    mat[i] = rr[i] * std::exp(c_i * TideFac::twopi() * uu[i]);
  }

  std::vector<int> ind = Constituent::iconst_unique();
  std::vector<std::complex<double>> F(Constituent::constituents()->size());
  std::vector<double> U(Constituent::constituents()->size());
  std::vector<double> V(Constituent::constituents()->size());
  for (size_t i = 0; i < ind.size(); ++i) {
    F[ind[i] - 1] = 1.0 + matsum(mat, ind[i]);
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

  std::vector<double> Fout(F.size());
  std::transform(F.begin(), F.end(), Fout.begin(),
                 [](std::complex<double> &d) { return d.real(); });

  return {Fout, U, V};
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

Date TideFac::curTime() const { return this->m_curTime; }

std::string TideFac::name(size_t index) {
  assert(index < this->m_tides.size());
  return this->m_tides[index].name;
}

double TideFac::amplitude(size_t index) {
  assert(index < this->m_tides.size());
  return this->m_tides[index].amp;
}

double TideFac::frequency(size_t index) {
  assert(index < this->m_tides.size());
  return this->m_tides[index].freq;
}

double TideFac::earthTideReductionFactor(size_t index) {
  assert(index < this->m_tides.size());
  return this->m_tides[index].etrf;
}

double TideFac::nodeFactor(size_t index) {
  assert(index < this->m_tides.size());
  return this->m_tides[index].nodefactor;
}

double TideFac::equilibriumArgument(size_t index) {
  assert(index < this->m_tides.size());
  return this->m_tides[index].eqarg;
}

Date TideFac::refTime() const { return this->m_refTime; }

void TideFac::setRefTime(const Date &refTime) { this->m_refTime = refTime; }
