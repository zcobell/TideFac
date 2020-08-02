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
#ifndef CONSTITUENT_H
#define CONSTITUENT_H

#include <array>
#include <string>
#include <limits>
#include <vector>

class Constituent {
 public:
  struct TC {
    const char *name;
    const double frequency;
    const std::array<double, 6> doodson;
    const double semi;
    const double isat;
    const int ishallow;
    const int nshallow;
    const double doodsonamp;
    const double earthreduc;
  };

  template <typename T>
  static constexpr T nullvalue() {
    return -std::numeric_limits<T>::max();
  }

  static size_t constituentIndex(const std::string &name);
  static const std::array<std::array<int, 3>, 162> *deldood();
  static const std::array<double, 162> *phcorr();
  static const std::array<double, 162> *amprat();
  static const std::array<int, 162> *ilatfac();
  static const std::array<int, 162> *iconst();
  static const std::array<int, 37> *iconst_unique();
  static const std::array<const TC *, 146> *constituents();
  static const std::array<int, 251> *shallow_iconst();
  static const std::array<double, 251> *shallow_coef();
  static const std::array<int, 251> *shallow_iname();
};

#endif  // CONSTITUENT_H
