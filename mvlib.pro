TEMPLATE = subdirs
CONFIG += ordered
CONFIG -= qt

SUBDIRS = mv \
          test

test.depends = mv
