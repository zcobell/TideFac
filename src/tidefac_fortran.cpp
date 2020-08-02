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
#include <vector>

#include "tidefac.h"

static std::vector<std::unique_ptr<TideFac>> s_tidefac;

extern "C" {
void* createTidefac();
void deleteTidefac(void* object);
int addConstituent(void* object, const char* name);
void setReferenceTime(void* object, int year, int month, int day, int hour,
                      int minute, int second);
void referenceTime(void* object, int& year, int& month, int& day, int& hour,
                   int& minute, int& second);
void calculateWithDt(void* object, double dt, const double latitude);
void calculateWithDate(void* object, int year, int month, int day, int hour,
                       int minute, int second, double latitude);
void calculateWithTwoDates(void* object, int year1, int month1, int day1,
                           int hour1, int minute1, int second1, int year2,
                           int month2, int day2, int hour2, int minute2,
                           int second2, double latitude);
void calculateGrid(void* object, double dt);
double amplitude(void* object, int index);
double frequency(void* object, int index);
double earthTideReductionFactor(void* object, int index);
double nodeFactor(void* object, int index);
double equilibriumArgument(void* object, int index);
double nodefactorCorrection(void* object, int index);
double astronomicArgument(void* object, int index);
double amplitudeGrid(void* object, int gridIndex, int index);
double frequencyGrid(void* object, int gridIndex, int index);
double earthTideReductionFactorGrid(void* object, int gridIndex, int index);
double nodeFactorGrid(void* object, int gridIndex, int index);
double equilibriumArgumentGrid(void* object, int gridIndex, int index);
double nodefactorCorrectionGrid(void* object, int gridIndex, int index);
double astronomicArgumentGrid(void* object, int gridIndex, int index);
void purgeTidefac();
int generateLatitudeGrid(void* object, double latmin, double latmax,
                         double resolution);
int getInterpolationFactors(void* object, double latitude,
                            unsigned long& gridIndex, double& weight);
}

void* createTidefac() {
  s_tidefac.push_back(std::unique_ptr<TideFac>(new TideFac()));
  return reinterpret_cast<void*>(s_tidefac.back().get());
}

int generateLatitudeGrid(void* object, double latmin, double latmax,
                         double resolution) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->generateLatitudeGrid(latmin, latmax, resolution);
}

int getInterpolationFactors(void* object, double latitude,
                            unsigned long& gridIndex, double& weight) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->getInterpolationFactors(latitude, gridIndex, weight);
}

void purgeTidefac() { s_tidefac.clear(); }

int addConstituent(void* object, const char* name) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->addConstituent(name);
}

void setReferenceTime(void* object, int year, int month, int day, int hour,
                      int minute, int second) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  f->setRefTime(Date(year, month, day, hour, minute, second));
}

void referenceTime(void* object, int& year, int& month, int& day, int& hour,
                   int& minute, int& second) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  Date d = f->refTime();
  year = d.year();
  month = d.month();
  day = d.day();
  hour = d.hour();
  minute = d.minute();
  second = d.second();
}

void calculateGrid(void* object, double dt) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  f->calculate(dt);
}

void calculateWithDt(void* object, double dt, double latitude) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  f->calculate(dt, latitude);
}

void calculateWithDate(void* object, int year, int month, int day, int hour,
                       int minute, int second, double latitude) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  f->calculate(Date(year, month, day, hour, minute, second), latitude);
}

void calculateWithTwoDates(void* object, int year1, int month1, int day1,
                           int hour1, int minute1, int second1, int year2,
                           int month2, int day2, int hour2, int minute2,
                           int second2, double latitude) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  f->calculate(Date(year1, month1, day1, hour1, minute1, second1),
               Date(year2, month2, day2, hour2, minute2, second2), latitude);
}

double amplitude(void* object, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->amplitude(index - 1);
}

double frequency(void* object, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->frequency(index - 1);
}

double earthTideReductionFactor(void* object, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->earthTideReductionFactor(index - 1);
}

double nodeFactor(void* object, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->nodeFactor(index - 1);
}

double equilibriumArgument(void* object, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->equilibriumArgument(index - 1);
}

double nodefactorCorrection(void* object, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->nodefactorCorrection(index - 1);
}

double astronomicArgument(void* object, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->astronomicArgument(index - 1);
}

double amplitudeGrid(void* object, int gridIndex, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->amplitude(index - 1, gridIndex - 1);
}

double frequencyGrid(void* object, int gridIndex, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->frequency(index - 1, gridIndex - 1);
}

double earthTideReductionFactorGrid(void* object, int gridIndex, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->earthTideReductionFactor(index - 1, gridIndex - 1);
}

double nodeFactorGrid(void* object, int gridIndex, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->nodeFactor(index - 1, gridIndex - 1);
}

double equilibriumArgumentGrid(void* object, int gridIndex, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->equilibriumArgument(index - 1, gridIndex - 1);
}

double nodefactorCorrectionGrid(void* object, int gridIndex, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->nodefactorCorrection(index - 1, gridIndex - 1);
}

double astronomicArgumentGrid(void* object, int gridIndex, int index) {
  TideFac* f = reinterpret_cast<TideFac*>(object);
  return f->astronomicArgument(index - 1, gridIndex - 1);
}
