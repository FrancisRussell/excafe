CXX=${MPICXX}
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc

bin_PROGRAMS = main solver

main_SOURCES = main.cpp
main_LDADD=$(top_srcdir)/lib/libsimplecfd.la

solver_SOURCES = solver.cpp
solver_LDADD=$(top_srcdir)/lib/libsimplecfd.la
