TEMPLATE = app
SOURCES = main.c
HEADERS = data.h
CONFIG -= qt
CONFIG += console

macx {
   CONFIG -= app_bundle
}

INCLUDEPATH = ../mv/include

LIBS += ../mv/libmv.a
TARGET = test_app

