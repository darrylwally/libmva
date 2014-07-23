# Filename: Makefile
# Description: Master makefile to make and build libmva.a and mvlib_test
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

all: libmva test 

.PHONY: libmva test libmva_debug test_debug

debug: libmva_debug test_debug

libmva:
	@cd mv && $(MAKE)

test: libmva
	@cd test && $(MAKE)

libmva_debug:
	@cd mv && $(MAKE) debug

test_debug: libmva_debug
	@cd test && $(MAKE) debug

clean: libmva_clean test_clean

libmva_clean:
	@cd mv && $(MAKE) clean

test_clean: 
	@cd test && $(MAKE) clean

distclean: libmva_distclean test_distclean

libmva_distclean:
	@cd mv && $(MAKE) distclean

test_distclean:
	@cd test && $(MAKE) distclean

first: all
