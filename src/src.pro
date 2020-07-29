#------------------------------GPL---------------------------------------//
# This file is part of TideFac.
#
# (c) 2020 Zachary Cobell
#
# TideFac is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TideFac is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with TideFac.  If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------*/
CONFIG -= qt

TEMPLATE = lib
DEFINES += TIDEFAC_LIBRARY
TARGET = tidefac

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    constituent.cpp \
    tidefac.cpp \
    date.cpp \
    tidefac_fortran.cpp

HEADERS += \
    TideFac_global.h \
    constituent.h \
    tidefac.h \
    date.h \
    date_hh.h

# Default rules for deployment.
unix {
    target.path = /usr/lib
}
!isEmpty(target.path): INSTALLS += target
