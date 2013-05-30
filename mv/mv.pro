TEMPLATE = lib
CONFIG += staticlib
CONFIG -= qt

include(cephes/cephes.pri)

INCLUDEPATH += include \
               ./
HEADERS += \
    include/mvstats.h \
    include/mvmodel.h \
    include/mvmatrix.h \
    include/mvconstants.h

SOURCES += \
    src/mvstats.c \
    src/mvmodel.c \
    src/mvmatrix.c \
    src/mvconstants.c

