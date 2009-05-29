CXX=${MPICXX}
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc

bin_PROGRAMS = tests
tests_SOURCES = tests.cpp tester.cpp
tests_LDADD=$(top_srcdir)/lib/libsimplecfd.la
