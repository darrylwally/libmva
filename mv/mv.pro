TEMPLATE = lib
CONFIG += staticlib
CONFIG -= qt

include(cephes/cephes.pri)

TARGET = mv
win32 {
    CONFIG(debug, debug|release) {
         TARGET = $$join(TARGET,,,d)
     }
    else{
        TARGET = $$TARGET
    }
}

DESTDIR = ./

INCLUDEPATH += include \
               ./
HEADERS += \
    include/mvstats.h \
    include/mvmodel.h \
    include/mvmatrix.h \
    include/mvconstants.h \
    include/mvpreprocess.h

SOURCES += \
    src/mvstats.c \
    src/mvmodel.c \
    src/mvmatrix.c \
    src/mvconstants.c \
    src/mvpreprocess.c

