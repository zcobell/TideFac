#include "tidefac.h"

#include <cmath>
#include <numeric>

TideFac::TideFac() {}

int TideFac::addConstituent(const char *harmonic) {
  return this->addConstituent(std::string(harmonic));
}

int TideFac::addConstituent(const std::string &harmonic) {
  auto c = Constituent(harmonic);
  if (c.isnull()) {
    return 1;
  } else {
    this->m_constituents.push_back(c);
    return 0;
  }
}

void TideFac::calculate(const Date &d, const double latitude) {
  this->computeAstronomicalArguments(d, latitude);
}

void TideFac::calculate(const size_t dt, const double latitude) {
  this->calculate(this->m_refTime + dt, latitude);
}

std::complex<double> matsum(const std::array<std::complex<double>, 162> &mat,
                            const int idx) {
  std::complex<double> s = {0.0, 0.0};
  for (size_t i = 0; i < 162; ++i) {
    if (Constituent::iconst()->at(i) == idx) {
      s += mat[i];
    }
  }
  return s;
}

std::tuple<std::vector<std::complex<double>>, std::vector<double>,
           std::vector<double>>
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

  const std::complex<double> c_i(0, 1);
  std::array<std::complex<double>, 162> mat;
  for (size_t i = 0; i < rr.size(); ++i) {
    mat[i] = rr[i] * std::exp(c_i * TideFac::twopi() * uu[i]);
  }

  std::vector<int> ind = Constituent::iconst_unique();
  std::vector<std::complex<double>> F(ind.size(), 1.0);
  std::vector<double> U(ind.size());
  std::vector<double> V(ind.size());
  for (size_t i = 0; i < ind.size(); ++i) {
    F[ind[i] - 1] = 1.0 + matsum(mat, ind[i]);
    U[i] = (std::log(F[i]).imag()) / (2.0 * TideFac::pi());
    F[i] = std::abs(F[i]);
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

      for (size_t p = 0; p < ik.size(); ++p) {
        double product_1 = 1.0;
        double product_2 = 0.0;
        for (size_t ii = 0; ii < ik.size(); ++ii) {
          product_1 = product_1 * std::pow(F[iloc[p] - 1].real(), exp2[ii]);
          product_2 += product_2 + U[iloc[p]] * exp1[ii];
        }
        F[i] = product_1;
        U[i] = product_2;
      }
    }
  }

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
