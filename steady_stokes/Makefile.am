CXX=mpicxx
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc -Wno-long-long

bin_PROGRAMS = main
main_SOURCES = main.cpp
main_LDADD=$(top_srcdir)/lib/libsimplecfd.la
