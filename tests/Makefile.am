CXX=${MPICXX}
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc ${PYTHON_CPPFLAGS}

bin_PROGRAMS = tests polynomials legendre expression_conversion sympy #forms

tests_SOURCES = tests.cpp tester.cpp
tests_LDADD=$(top_srcdir)/lib/libsimplecfd.la

polynomials_SOURCES = polynomials.cpp
polynomials_LDADD=$(top_srcdir)/lib/libsimplecfd.la

legendre_SOURCES = legendre.cpp
legendre_LDADD=$(top_srcdir)/lib/libsimplecfd.la

expression_conversion_SOURCES = expression_conversion.cpp
expression_conversion_LDADD=$(top_srcdir)/lib/libsimplecfd.la

sympy_SOURCES = sympy.cpp
sympy_LDADD=$(top_srcdir)/lib/libsimplecfd.la

#forms_SOURCES = forms.cpp
#forms_LDADD=$(top_srcdir)/lib/libsimplecfd.la
