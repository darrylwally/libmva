# Filename: Makefile
# Description: Master makefile to make and build libmv.a and mvlib_test

all: libmv test

.PHONY: libmv test libmv_debug test_debug

debug: libmv_debug test_debug

libmv:
	@cd mv && $(MAKE)

test:
	@cd test && $(MAKE)

libmv_debug:
	@cd mv && $(MAKE) debug

test_debug:
	@cd test && $(MAKE) debug

clean: libmv_clean test_clean

libmv_clean:
	@cd mv && $(MAKE) clean

test_clean: 
	@cd test && $(MAKE) clean

distclean: libmv_distclean test_distclean

libmv_distclean:
	@cd mv && $(MAKE) distclean

test_distclean:
	@cd test && $(MAKE) distclean

first: all
