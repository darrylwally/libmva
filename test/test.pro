TEMPLATE = app
SOURCES = main.c
HEADERS = data.h
CONFIG -= qt
CONFIG += console

macx {
   CONFIG -= app_bundle
}

INCLUDEPATH = ../mv/include

win32: {
    LIBS += Psapi.lib
    DEFINES +=PSAPI_VERSION=1
    QMAKE_CFLAGS += -WX
    CONFIG(debug, debug|release) {
         LIBS += ../mv/mvd.lib
     }
    else{
        LIBS += ../mv/mv.lib
    }
}
else: {
LIBS += ../mv/libmv.a
}
TARGET = test_app

