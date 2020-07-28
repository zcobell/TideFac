/*------------------------------GPL---------------------------------------//
// This file is part of HMDF.
//
// (c) 2015-2020 Zachary Cobell
//
// HMDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HMDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with HMDF.  If not, see <http://www.gnu.org/licenses/>.
//------------------------------------------------------------------------*/
#include "date.h"

#include <cassert>
#include <chrono>
#include <iostream>

#include "date_hh.h"

struct s_date {
 private:
  date::year_month_day dd;

 public:
  explicit s_date(const std::chrono::system_clock::time_point &t)
      : dd(date::year_month_day(date::floor<date::days>(t))) {}
  int year() const { return int(dd.year()); }
  unsigned month() const { return unsigned(dd.month()); }
  unsigned day() const { return unsigned(dd.day()); }
};

struct s_datetime {
 private:
  std::chrono::system_clock::time_point t;
  date::year_month_day dd;
  date::hh_mm_ss<std::chrono::system_clock::duration> tt;

 public:
  explicit s_datetime(const std::chrono::system_clock::time_point &t)
      : t(t),
        dd(date::year_month_day(date::floor<date::days>(t))),
        tt(date::make_time(t - date::sys_days(dd))) {}
  date::year_month_day ymd() { return dd; }
  int year() const { return int(dd.year()); }
  unsigned month() const { return unsigned(dd.month()); }
  unsigned day() const { return unsigned(dd.day()); }
  int hour() const { return tt.hours().count(); }
  int minute() const { return tt.minutes().count(); }
  int second() const { return tt.seconds().count(); }
  int milliseconds() const {
    Date c(year(), month(), day(), hour(), minute(), second());
    return std::chrono::duration_cast<std::chrono::milliseconds>(t -
                                                                 c.time_point())
        .count();
  }
};

date::year_month_day normalize(date::year_month_day ymd) {
  ymd += date::months{0};
  ymd = date::sys_days{ymd};
  return ymd;
}

constexpr date::year_month_day c_epoch() {
  return date::year_month_day(date::year(1970) / 1 / 1);
}

Date::Date() { this->set(1970, 1, 1, 0, 0, 0, 0); }

Date::Date(const std::chrono::system_clock::time_point &t) { this->set(t); }

Date::Date(const std::vector<int> &v) { this->set(v); }

Date::Date(const Date &d) { this->set(d.get()); }

Date::Date(int year, int month, int day, int hour, int minute, int second,
           int millisecond) {
  this->set(year, month, day, hour, minute, second, millisecond);
}

void Date::addSeconds(const long &value) { this->m_datetime += seconds(value); }

void Date::addMinutes(const long &value) { this->m_datetime += minutes(value); }

void Date::addHours(const long &value) { this->m_datetime += hours(value); }

void Date::addDays(const long &value) { this->m_datetime += days(value); }

void Date::addWeeks(const long &value) { this->m_datetime += weeks(value); }

void Date::addMonths(const long &value) { this->m_datetime += months(value); }

void Date::addYears(const long &value) { this->m_datetime += years(value); }

bool Date::isLeap() const { return date::year(this->year()).is_leap(); }

bool Date::operator<(const Date &d) const {
  return this->time_point() < d.time_point();
}

bool Date::operator>(const Date &d) const {
  return this->time_point() > d.time_point();
}

bool Date::operator==(const Date &d) const {
  return this->time_point() == d.time_point();
}

bool Date::operator!=(const Date &d) const { return !(*(this) == d); }

std::ostream &operator<<(std::ostream &os, const Date &dt) {
  os << dt.toString();
  return os;
}

Date &Date::operator-=(const Date::months &rhs) {
  s_datetime d(this->m_datetime);
  date::year_month_day dd = d.ymd();
  dd -= date::months(rhs);
  dd = normalize(dd);
  this->set(int(dd.year()), unsigned(dd.month()), unsigned(dd.day()),
            this->hour(), this->minute(), this->second(), this->millisecond());
  return *this;
}

Date &Date::operator-=(const Date::years &rhs) {
  s_datetime d(this->m_datetime);
  date::year_month_day dd = d.ymd();
  dd -= date::years(rhs);
  dd = normalize(dd);
  this->set(int(dd.year()), unsigned(dd.month()), unsigned(dd.day()),
            this->hour(), this->minute(), this->second(), this->millisecond());
  return *this;
}

Date &Date::operator+=(const Date::months &rhs) {
  s_datetime d(this->m_datetime);
  date::year_month_day dd = d.ymd();
  dd += date::months(rhs);
  dd = normalize(dd);
  this->set(int(dd.year()), unsigned(dd.month()), unsigned(dd.day()),
            this->hour(), this->minute(), this->second(), this->millisecond());
  return *this;
}

Date &Date::operator+=(const Date::years &rhs) {
  s_datetime d(this->m_datetime);
  date::year_month_day dd = d.ymd();
  dd += date::years(rhs);
  dd = normalize(dd);
  this->set(int(dd.year()), unsigned(dd.month()), unsigned(dd.day()),
            this->hour(), this->minute(), this->second(), this->millisecond());
  return *this;
}

void Date::set(int year, int month, int day, int hour, int minute, int second,
               int millisecond) {
  auto ymd = date::year(year) / date::month(month) / date::day(day);
  if (!ymd.ok()) {
    // hmdf_throw_exception("Invalid date");
    return;
  }
  this->m_datetime = date::sys_days(ymd) + std::chrono::hours(hour) +
                     std::chrono::minutes(minute) +
                     std::chrono::seconds(second) +
                     std::chrono::milliseconds(millisecond);
}

std::vector<int> Date::get() const {
  s_datetime time(this->m_datetime);
  std::vector<int> v(7);
  v[0] = time.year();
  v[1] = time.month();
  v[2] = time.day();
  v[3] = time.hour();
  v[4] = time.minute();
  v[5] = time.second();
  v[6] = time.milliseconds();
  return v;
}

void Date::set(const std::vector<int> &v) {
  std::vector<int> v2(7);
  std::copy(v.begin(), v.end(), v2.begin());

  this->set(v[0], v[1], v[2], v[3], v[4], v[5], v[6]);
}

void Date::set(const std::chrono::system_clock::time_point &t) {
  this->m_datetime = t;
}

void Date::set(const Date &v) { this->set(v.time_point()); }

void Date::fromSeconds(long seconds) {
  this->m_datetime = date::sys_days(c_epoch()) + std::chrono::seconds(seconds);
}

void Date::fromMSeconds(long long mseconds) {
  this->fromSeconds(mseconds / 1000);
}

long Date::toSeconds() const {
  return std::chrono::duration_cast<std::chrono::seconds>(
             this->m_datetime - date::sys_days(c_epoch()))
      .count();
}

long long Date::toMSeconds() const { return this->toSeconds() * 1000; }

int Date::year() const {
  s_date d(this->m_datetime);
  return d.year();
}

void Date::setYear(int year) {
  s_datetime d(this->m_datetime);
  this->set(year, d.month(), d.day(), d.hour(), d.minute(), d.second());
}

int Date::month() const {
  s_date d(this->m_datetime);
  return d.month();
}

void Date::setMonth(int month) {
  s_datetime d(this->m_datetime);
  this->set(d.year(), month, d.day(), d.hour(), d.minute(), d.second());
}

int Date::day() const {
  s_date d(this->m_datetime);
  return d.day();
}

void Date::setDay(int day) {
  s_datetime d(this->m_datetime);
  this->set(d.year(), d.month(), day, d.hour(), d.minute(), d.second());
}

int Date::hour() const {
  s_datetime d(this->m_datetime);
  return d.hour();
}

void Date::setHour(int hour) {
  s_datetime d(this->m_datetime);
  this->set(d.year(), d.month(), d.day(), hour, d.minute(), d.second());
}

int Date::minute() const {
  s_datetime d(this->m_datetime);
  return d.minute();
}

void Date::setMinute(int minute) {
  s_datetime d(this->m_datetime);
  this->set(d.year(), d.month(), d.day(), d.hour(), minute, d.second());
}

int Date::second() const {
  s_datetime d(this->m_datetime);
  return d.second();
}

void Date::setSecond(int second) {
  s_datetime d(this->m_datetime);
  this->set(d.year(), d.month(), d.day(), d.hour(), d.minute(), second);
}

int Date::millisecond() const {
  s_datetime d(this->m_datetime);
  return d.milliseconds();
}

void Date::setMillisecond(int milliseconds) {
  s_datetime d(this->m_datetime);
  this->set(d.year(), d.month(), d.day(), d.hour(), d.minute(), d.second(),
            milliseconds);
}

void Date::fromString(const std::string &datestr, const std::string &format) {
  std::stringstream ss(datestr);
  date::from_stream(ss, format.c_str(), this->m_datetime);
}

std::string Date::toString(const std::string &format) const {
  return date::format(format, this->m_datetime);
}

std::chrono::system_clock::time_point Date::time_point() const {
  return this->m_datetime;
}

Date Date::now() { return Date(std::chrono::system_clock::now()); }
