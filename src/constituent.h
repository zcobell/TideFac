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

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "constituent_data.h"

namespace Constituent {

constexpr size_t constituentIndex(const std::string_view &name) {
  for (size_t i = 0; i < s_names.size(); ++i) {
    if (s_names[i] == name) {
      return i;
    }
  }
  return Constituent::null_value<size_t>();
}

constexpr std::array<std::array<char, 3>, 162> deldood() { return s_deldood; }

constexpr std::array<double, 162> phcorr() { return s_phcorr; }

constexpr std::array<double, 162> amprat() { return s_amprat; }

constexpr std::array<unsigned char, 162> ilatfac() { return s_ilatfac; }

constexpr std::array<unsigned char, 162> iconst() { return s_iconst; }

constexpr std::array<unsigned char, 251> shallow_iconst() {
  return s_shallow_iconst;
}

constexpr std::array<double, 251> shallow_coef() { return s_shallow_coef; }

constexpr std::array<unsigned char, 251> shallow_iname() {
  return s_shallow_iname;
}

constexpr std::array<Constituent::TC, 146> constituents() { return s_alltides; }

constexpr std::array<short, 37> iconst_unique() { return s_iconst_unique; }
};  // namespace Constituent

#endif  // CONSTITUENT_H
