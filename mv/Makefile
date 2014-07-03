# Filename: Makefile
# Description: Makefile to build static library for libmv
# Author: Darryl Wallace <wallacdj@gmail.com
# Copyright (c) 2014 - Darryl Wallace
#
# License:
# libmva is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# libmva is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

##### Compiler tools, etc.
AR            = ar cq
RANLIB        = ranlib -s
CC            = clang
CFLAGS        = -pipe -gdwarf-2 -fPIC -arch x86_64 -Wall
INCPATH       = -I. -Iinclude
DEL_FILE      = rm -f
DEL_DIR       = rmdir
COPY_FILE     = cp -f

##### Files
SOURCES       = cephes/cmath/const.c \
                cephes/cmath/isnan.c \
                cephes/cmath/mtherr.c \
                cephes/cprob/chdtr.c \
                cephes/cprob/fdtr.c \
                cephes/cprob/gamma.c \
                cephes/cprob/incbet.c \
                cephes/cprob/incbi.c \
                cephes/cprob/igam.c \
                cephes/cprob/igami.c \
                cephes/cprob/ndtr.c \
                cephes/cprob/ndtri.c \
                cephes/cprob/polevl.c \
                cephes/cprob/stdtr.c \
                cephes/misc/fsolve.c \
                cephes/misc/gammaincinv.c \
                src/mvstats.c \
                src/mvmodel.c \
                src/mvmatrix.c \
                src/mvconstants.c \
                src/mvpreprocess.c \
                src/mvfile.c

OBJECTS       = const.o \
                isnan.o \
                mtherr.o \
                chdtr.o \
                fdtr.o \
                gamma.o \
                incbet.o \
                incbi.o \
                igam.o \
                igami.o \
                ndtr.o \
                ndtri.o \
                polevl.o \
                stdtr.o \
                fsolve.o \
                gammaincinv.o \
                mvstats.o \
                mvmodel.o \
                mvmatrix.o \
                mvconstants.o \
                mvpreprocess.o \
                mvfile.o

TARGET        = libmva.a
DEBUG_TARGET  = libmvad.a

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all:    release

#staticlib: $(TARGET)

release: CFLAGS += -O3
release: staticlib

debug: CFLAGS += -g
debug: TARGET = $(DEBUG_TARGET)
debug: staticlib

staticlib:  $(OBJECTS) $(OBJCOMP)
	-$(DEL_FILE) $(TARGET)	
	$(AR) $(TARGET) $(OBJECTS)
	$(RANLIB) $(TARGET)

####### Compile

const.o: cephes/cmath/const.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o const.o cephes/cmath/const.c

isnan.o: cephes/cmath/isnan.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o isnan.o cephes/cmath/isnan.c

mtherr.o: cephes/cmath/mtherr.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o mtherr.o cephes/cmath/mtherr.c

chdtr.o: cephes/cprob/chdtr.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o chdtr.o cephes/cprob/chdtr.c

fdtr.o: cephes/cprob/fdtr.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o fdtr.o cephes/cprob/fdtr.c

gamma.o: cephes/cprob/gamma.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o gamma.o cephes/cprob/gamma.c

incbet.o: cephes/cprob/incbet.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o incbet.o cephes/cprob/incbet.c

incbi.o: cephes/cprob/incbi.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o incbi.o cephes/cprob/incbi.c

igam.o: cephes/cprob/igam.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o igam.o cephes/cprob/igam.c

igami.o: cephes/cprob/igami.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o igami.o cephes/cprob/igami.c

ndtr.o: cephes/cprob/ndtr.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o ndtr.o cephes/cprob/ndtr.c

ndtri.o: cephes/cprob/ndtri.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o ndtri.o cephes/cprob/ndtri.c

polevl.o: cephes/cprob/polevl.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o polevl.o cephes/cprob/polevl.c

stdtr.o: cephes/cprob/stdtr.c cephes/cmath/mconf.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o stdtr.o cephes/cprob/stdtr.c

fsolve.o: cephes/misc/fsolve.c cephes/misc/misc.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o fsolve.o cephes/misc/fsolve.c

gammaincinv.o: cephes/misc/gammaincinv.c cephes/cmath/mconf.h \
                cephes/misc/misc.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o gammaincinv.o cephes/misc/gammaincinv.c

mvstats.o: src/mvstats.c cephes/cmath/mconf.h \
                cephes/cmath/protos.h
	$(CC) -c $(CFLAGS) $(INCPATH) -o mvstats.o src/mvstats.c

mvmodel.o: src/mvmodel.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o mvmodel.o src/mvmodel.c

mvconstants.o: src/mvconstants.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o mvconstants.o src/mvconstants.c

mvpreprocess.o: src/mvpreprocess.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o mvpreprocess.o src/mvpreprocess.c

mvmatrix.o: src/mvmatrix.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o mvmatrix.o src/mvmatrix.c

mvfile.o: src/mvmatrix.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o mvfile.o src/mvfile.c

##### Clean
clean: compiler_clean
	-$(DEL_FILE) $(OBJECTS)

distclean: clean
	-$(DEL_FILE) $(TARGET)
	-$(DEL_FILE) $(DEBUG_TARGET)

compiler_clean:

##### Misc
check: first
