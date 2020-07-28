QT -= gui
QT -= qt

CONFIG += c++11 console
CONFIG -= app_bundle

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
        bench.cpp

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../src/release/ -ltidefac
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../src/debug/ -ltidefac
else:unix: LIBS += -L$$OUT_PWD/../src/ -ltidefac

INCLUDEPATH += $$PWD/../src
DEPENDPATH += $$PWD/../src


macx: LIBS += -L/opt/google-benchmark/lib/ -lbenchmark

INCLUDEPATH += /opt/google-benchmark/include
DEPENDPATH += /opt/google-benchmark/include

macx: PRE_TARGETDEPS += /opt/google-benchmark/lib/libbenchmark.a
