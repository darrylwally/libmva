TEMPLATE = app
CONFIG += console
CONFIG -= qt


SOURCES += \
    main.c \
    mvmodel.c \
    mvconstants.c \
    mvmatrix.c \
    mvstats.c \
    cephes/besselpoly.c \
    cephes/chdtr.c \
    cephes/const.c \
    cephes/fdtr.c \
    cephes/fsolve.c \
    cephes/gamma.c \
    cephes/gammaincinv.c \
    cephes/igami.c \
    cephes/igam.c \
    cephes/incbet.c \
    cephes/incbi.c \
    cephes/mtherr.c \
    cephes/ndtri.c \
    cephes/polevl.c \
    cephes/stdtr.c

HEADERS += \
    data.h \
    mvmodel.h \
    mvconstants.h \
    mvmatrix.h \
    mvstats.h \
    cephes/mconf.h \
    cephes/protos.h

win32: {
    LIBS += Psapi.lib
    DEFINES +=PSAPI_VERSION=1
    QMAKE_CFLAGS += -WX
}
unix: {
    QMAKE_CFLAGS += -Werror
}

