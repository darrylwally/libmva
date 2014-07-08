# Filename: Makefile
# Description: Makefile to test program for libmv
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
CLFAGS        = -pipe -gdwarf-2 -fPIC -arch x86_64 -Wall
INCPATH       = -I. -I../mv/include
LIBS          = ../mv/libmv.a
LINK          = clang++
LFLAGS        = -arch x86_64
DEL_FILE      = rm -f
DEL_DIR       = rmdir
COPY_FILE     = cp -f

##### Files
SOURCES       = main.c
OBJECTS       = main.o

TARGET        = mvlib_test
DEBUG_TARGET  = mvlib_testd

first: all
##### Implicit rules
.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

##### Build Rules

all: release

release: CFLAGS += -O3
release: executable

debug: CFLAGS += -g
debug: LIBS = ../mv/libmvd.a
debug: TARGET = $(DEBUG_TARGET)
debug: executable

executable: $(OBJECTS)
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS)

##### Clean
clean: compiler_clean
	-$(DEL_FILE) $(OBJECTS)

distclean: clean
	-$(DEL_FILE) $(TARGET)
	-$(DEL_FILE) $(DEBUG_TARGET)

compiler_clean:

##### Compile

main.o: main.c data.h
	$(CC) -c $(CLFAGS) $(INCPATH) -o main.o main.c
