CXX=${MPICXX}
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc

bin_PROGRAMS = cse tests polynomials legendre expression_conversion #forms

cse_SOURCES = cse.cpp
cse_LDADD=$(top_srcdir)/lib/libsimplecfd.la

tests_SOURCES = tests.cpp tester.cpp
tests_LDADD=$(top_srcdir)/lib/libsimplecfd.la

polynomials_SOURCES = polynomials.cpp
polynomials_LDADD=$(top_srcdir)/lib/libsimplecfd.la

legendre_SOURCES = legendre.cpp
legendre_LDADD=$(top_srcdir)/lib/libsimplecfd.la

expression_conversion_SOURCES = expression_conversion.cpp
expression_conversion_LDADD=$(top_srcdir)/lib/libsimplecfd.la

#forms_SOURCES = forms.cpp
#forms_LDADD=$(top_srcdir)/lib/libsimplecfd.la
