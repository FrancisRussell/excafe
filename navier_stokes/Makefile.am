CXX=${MPICXX}
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc

bin_PROGRAMS = solver

solver_SOURCES = solver.cpp
solver_LDADD=$(top_srcdir)/lib/libsimplecfd.la
