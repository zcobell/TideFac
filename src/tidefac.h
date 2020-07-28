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

  void TIDEFAC_EXPORT show() const;

  Date TIDEFAC_EXPORT curTime() const;

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
