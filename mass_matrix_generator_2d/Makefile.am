CXX=${MPICXX}
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc

bin_PROGRAMS = generator

generator_SOURCES = generator.cpp
generator_LDADD=$(top_srcdir)/lib/libsimplecfd.la
