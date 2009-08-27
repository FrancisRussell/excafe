CXX=${MPICXX}
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc

bin_PROGRAMS = tests polynomials legendre #forms

tests_SOURCES = tests.cpp tester.cpp
tests_LDADD=$(top_srcdir)/lib/libsimplecfd.la

polynomials_SOURCES = polynomials.cpp
polynomials_LDADD=$(top_srcdir)/lib/libsimplecfd.la

legendre_SOURCES = legendre.cpp
legendre_LDADD=$(top_srcdir)/lib/libsimplecfd.la

#forms_SOURCES = forms.cpp
#forms_LDADD=$(top_srcdir)/lib/libsimplecfd.la
