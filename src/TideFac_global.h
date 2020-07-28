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
#ifndef TIDEFAC_GLOBAL_H
#define TIDEFAC_GLOBAL_H

#if defined(_MSC_VER) || defined(WIN64) || defined(_WIN64) ||  \
    defined(__WIN64__) || defined(WIN32) || defined(_WIN32) || \
    defined(__WIN32__) || defined(__NT__)
#define Q_DECL_EXPORT __declspec(dllexport)
#define Q_DECL_IMPORT __declspec(dllimport)
#elif defined(SWIG)
#define Q_DECL_EXPORT
#define Q_DECL_IMPORT
#else
#define Q_DECL_EXPORT __attribute__((visibility("default")))
#define Q_DECL_IMPORT __attribute__((visibility("default")))
#endif

#if defined(TIDEFAC_LIBRARY)
#define TIDEFAC_EXPORT Q_DECL_EXPORT
#else
#define TIDEFAC_EXPORT Q_DECL_IMPORT
#endif

#endif  // TIDEFAC_GLOBAL_H
