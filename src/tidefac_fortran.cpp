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
#include <memory>
#include <unordered_map>

#include "tidefac.h"

static std::unordered_map<size_t, std::unique_ptr<TideFac>> s_tidefac;
static size_t tid = -1;

extern "C" {
size_t createTidefac();
void deleteTidefac(size_t object);
int addConstituent(size_t object, const char* name);
void setReferenceTime(size_t object, int year, int month, int day, int hour,
                      int minute, int second);
void referenceTime(size_t object, int& year, int& month, int& day, int& hour,
                   int& minute, int& second);
void calculateWithDt(size_t object, double dt, const double latitude);
void calculateWithDate(size_t object, int year, int month, int day, int hour,
                       int minute, int second, double latitude);
void calculateWithTwoDates(size_t object, int year1, int month1, int day1,
                           int hour1, int minute1, int second1, int year2,
                           int month2, int day2, int hour2, int minute2,
                           int second2, double latitude);
void calculateGrid(size_t object, double dt);
double amplitude(size_t object, int index);
double frequency(size_t object, int index);
double earthTideReductionFactor(size_t object, int index);
double nodeFactor(size_t object, int index);
double equilibriumArgument(size_t object, int index);
double nodefactorCorrection(size_t object, int index);
double astronomicArgument(size_t object, int index);
double amplitudeGrid(size_t object, int gridIndex, int index);
double frequencyGrid(size_t object, int gridIndex, int index);
double earthTideReductionFactorGrid(size_t object, int gridIndex, int index);
double nodeFactorGrid(size_t object, int gridIndex, int index);
double equilibriumArgumentGrid(size_t object, int gridIndex, int index);
double nodefactorCorrectionGrid(size_t object, int gridIndex, int index);
double astronomicArgumentGrid(size_t object, int gridIndex, int index);
void purgeTidefac();
int generateLatitudeGrid(size_t object, double latmin, double latmax,
                         double resolution);
int getInterpolationFactors(size_t object, double latitude,
                            unsigned long& gridIndex, double& weight);
}

size_t createTidefac() {
  tid++;
  s_tidefac[tid] = std::make_unique<TideFac>();
  return tid;
}

void deleteTidefac(size_t object) { s_tidefac.erase(object); }

int generateLatitudeGrid(size_t object, double latmin, double latmax,
                         double resolution) {
  return s_tidefac[object]->generateLatitudeGrid(latmin, latmax, resolution);
}

int getInterpolationFactors(size_t object, double latitude,
                            unsigned long& gridIndex, double& weight) {
  return s_tidefac[object]->getInterpolationFactors(latitude, gridIndex,
                                                    weight);
}

void purgeTidefac() { s_tidefac.clear(); }

int addConstituent(size_t object, const char* name) {
  return s_tidefac[object]->addConstituent(name);
}

void setReferenceTime(size_t object, int year, int month, int day, int hour,
                      int minute, int second) {
  s_tidefac[object]->setRefTime(Date(year, month, day, hour, minute, second));
}

void referenceTime(size_t object, int& year, int& month, int& day, int& hour,
                   int& minute, int& second) {
  Date d = s_tidefac[object]->refTime();
  year = d.year();
  month = d.month();
  day = d.day();
  hour = d.hour();
  minute = d.minute();
  second = d.second();
}

void calculateGrid(size_t object, double dt) {
  s_tidefac[object]->calculate(dt);
}

void calculateWithDt(size_t object, double dt, double latitude) {
  s_tidefac[object]->calculate(dt, latitude);
}

void calculateWithDate(size_t object, int year, int month, int day, int hour,
                       int minute, int second, double latitude) {
  s_tidefac[object]->calculate(Date(year, month, day, hour, minute, second),
                               latitude);
}

void calculateWithTwoDates(size_t object, int year1, int month1, int day1,
                           int hour1, int minute1, int second1, int year2,
                           int month2, int day2, int hour2, int minute2,
                           int second2, double latitude) {
  s_tidefac[object]->calculate(
      Date(year1, month1, day1, hour1, minute1, second1),
      Date(year2, month2, day2, hour2, minute2, second2), latitude);
}

double amplitude(size_t object, int index) {
  return s_tidefac[object]->amplitude(index - 1);
}

double frequency(size_t object, int index) {
  return s_tidefac[object]->frequency(index - 1);
}

double earthTideReductionFactor(size_t object, int index) {
  return s_tidefac[object]->earthTideReductionFactor(index - 1);
}

double nodeFactor(size_t object, int index) {
  return s_tidefac[object]->nodeFactor(index - 1);
}

double equilibriumArgument(size_t object, int index) {
  return s_tidefac[object]->equilibriumArgument(index - 1);
}

double nodefactorCorrection(size_t object, int index) {
  return s_tidefac[object]->nodefactorCorrection(index - 1);
}

double astronomicArgument(size_t object, int index) {
  return s_tidefac[object]->astronomicArgument(index - 1);
}

double amplitudeGrid(size_t object, int gridIndex, int index) {
  return s_tidefac[object]->amplitude(index - 1, gridIndex - 1);
}

double frequencyGrid(size_t object, int gridIndex, int index) {
  return s_tidefac[object]->frequency(index - 1, gridIndex - 1);
}

double earthTideReductionFactorGrid(size_t object, int gridIndex, int index) {
  return s_tidefac[object]->earthTideReductionFactor(index - 1, gridIndex - 1);
}

double nodeFactorGrid(size_t object, int gridIndex, int index) {
  return s_tidefac[object]->nodeFactor(index - 1, gridIndex - 1);
}

double equilibriumArgumentGrid(size_t object, int gridIndex, int index) {
  return s_tidefac[object]->equilibriumArgument(index - 1, gridIndex - 1);
}

double nodefactorCorrectionGrid(size_t object, int gridIndex, int index) {
  return s_tidefac[object]->nodefactorCorrection(index - 1, gridIndex - 1);
}

double astronomicArgumentGrid(size_t object, int gridIndex, int index) {
  return s_tidefac[object]->astronomicArgument(index - 1, gridIndex - 1);
}
