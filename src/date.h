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
#ifndef TDATE_H
#define TDATE_H

#include <chrono>
#include <cmath>
#include <iostream>
#include <ratio>
#include <string>
#include <type_traits>
#include <vector>

class Date {
 public:
  using milliseconds = std::chrono::milliseconds;
  using seconds = std::chrono::seconds;
  using minutes = std::chrono::minutes;
  using hours = std::chrono::hours;
  using days = std::chrono::duration<
      int, std::ratio_multiply<std::ratio<24>, std::chrono::hours::period>>;
  using weeks =
      std::chrono::duration<int,
                            std::ratio_multiply<std::ratio<7>, days::period>>;
  using years = std::chrono::duration<
      int, std::ratio_multiply<std::ratio<146097, 400>, days::period>>;
  using months =
      std::chrono::duration<int,
                            std::ratio_divide<years::period, std::ratio<12>>>;

  Date();
  Date(const std::chrono::system_clock::time_point &t);
  Date(const std::vector<int> &v);
  Date(const Date &d);
  Date(int year, int month = 1, int day = 1, int hour = 0, int minute = 0,
       int second = 0, int millisecond = 0);

#ifndef SWIG
  //...operator overloads
  bool operator<(const Date &d) const;
  bool operator>(const Date &d) const;
  bool operator==(const Date &d) const;
  bool operator!=(const Date &d) const;

  template <class T, typename std::enable_if<std::is_integral<T>::value>::type
                         * = nullptr>
  Date &operator+=(const T &rhs) {
    this->m_datetime += Date::seconds(rhs);
    return *this;
  }

  template <class T, typename std::enable_if<
                         std::is_floating_point<T>::value>::type * = nullptr>
  Date &operator+=(const T &rhs) {
    this->m_datetime +=
        Date::milliseconds(static_cast<long>(std::floor(rhs * 1000.0)));
    return *this;
  }

  template <class T,
            typename std::enable_if<!std::is_integral<T>::value &&
                                    !std::is_floating_point<T>::value &&
                                    !std::is_same<T, Date::years>::value &&
                                    !std::is_same<T, Date::months>::value>::type
                * = nullptr>
  Date &operator+=(const T &rhs) {
    this->m_datetime += rhs;
    return *this;
  }

  Date &operator+=(const Date::years &rhs);
  Date &operator+=(const Date::months &rhs);

  template <class T, typename std::enable_if<std::is_integral<T>::value>::type
                         * = nullptr>
  Date &operator-=(const T &rhs) {
    this->m_datetime -= Date::seconds(rhs);
    return *this;
  }

  template <class T, typename std::enable_if<
                         std::is_floating_point<T>::value>::type * = nullptr>
  Date &operator-=(const T &rhs) {
    this->m_datetime -=
        Date::milliseconds(static_cast<long>(std::floor(rhs * 1000.0)));
    return *this;
  }

  template <class T,
            typename std::enable_if<!std::is_integral<T>::value &&
                                    !std::is_floating_point<T>::value &&
                                    !std::is_same<T, Date::years>::value &&
                                    !std::is_same<T, Date::months>::value>::type
                * = nullptr>
  Date &operator-=(const T &rhs) {
    this->m_datetime -= rhs;
    return *this;
  }

  Date &operator-=(const Date::years &rhs);
  Date &operator-=(const Date::months &rhs);

  friend std::ostream &operator<<(std::ostream &os, const Date &dt);

  // end operator overloads
#endif

  void addSeconds(const long &value);
  void addMinutes(const long &value);
  void addHours(const long &value);
  void addDays(const long &value);
  void addWeeks(const long &value);
  void addMonths(const long &value);
  void addYears(const long &value);
  bool isLeap() const;

  static Date maxDate() { return Date(3000, 1, 1, 0, 0, 0); }
  static Date minDate() { return Date(1900, 1, 1, 0, 0, 0); }

  std::vector<int> get() const;

  void set(const std::vector<int> &v);
  void set(const std::chrono::system_clock::time_point &t);
  void set(const Date &v);
  void set(int year, int month = 1, int day = 1, int hour = 0, int minute = 0,
           int second = 0, int millisecond = 0);

  void fromSeconds(long seconds);

  void fromMSeconds(long long mseconds);

  long toSeconds() const;

  long long toMSeconds() const;

  int year() const;
  void setYear(int year);

  int month() const;
  void setMonth(int month);

  int day() const;
  void setDay(int day);

  int hour() const;
  void setHour(int hour);

  int minute() const;
  void setMinute(int minute);

  int second() const;
  void setSecond(int second);

  int millisecond() const;
  void setMillisecond(int milliseconds);

  void fromString(const std::string &datestr,
                  const std::string &format = "%Y-%m-%d %H:%M:%OS");

  std::string toString(const std::string &format = "%Y-%m-%d %H:%M:%OS") const;

  std::chrono::system_clock::time_point time_point() const;

  static Date now();

 private:
  std::chrono::system_clock::time_point m_datetime;
};

template <typename T>
Date operator+(Date lhs, const T &rhs) {
  lhs += rhs;
  return lhs;
}

template <typename T>
Date operator-(Date lhs, const T &rhs) {
  lhs -= rhs;
  return lhs;
}

#ifndef SWIG
template Date operator+(Date, const short &);
template Date operator+(Date, const int &);
template Date operator+(Date, const long &);
template Date operator+(Date, const unsigned short &);
template Date operator+(Date, const unsigned int &);
template Date operator+(Date, const unsigned long &);
template Date operator+(Date, const float &);
template Date operator+(Date, const double &);
template Date operator+(Date, const Date::milliseconds &);
template Date operator+(Date, const Date::seconds &);
template Date operator+(Date, const Date::minutes &);
template Date operator+(Date, const Date::hours &);
template Date operator+(Date, const Date::days &);
template Date operator+(Date, const Date::months &);
template Date operator+(Date, const Date::weeks &);
template Date operator+(Date, const Date::years &);
template Date operator-(Date, const short &);
template Date operator-(Date, const int &);
template Date operator-(Date, const long &);
template Date operator-(Date, const unsigned short &);
template Date operator-(Date, const unsigned int &);
template Date operator-(Date, const unsigned long &);
template Date operator-(Date, const float &);
template Date operator-(Date, const double &);
template Date operator-(Date, const Date::milliseconds &);
template Date operator-(Date, const Date::seconds &);
template Date operator-(Date, const Date::minutes &);
template Date operator-(Date, const Date::hours &);
template Date operator-(Date, const Date::days &);
template Date operator-(Date, const Date::months &);
template Date operator-(Date, const Date::weeks &);
template Date operator-(Date, const Date::years &);
#endif

#endif  // TDATE_H
