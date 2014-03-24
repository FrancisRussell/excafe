CXX=${MPICXX}
AM_CXXFLAGS=-I$(top_srcdir)/include -I/usr/include/petsc

bin_PROGRAMS = generator

generator_SOURCES = generator.cpp mat_mult_code_generator.cpp vector_entry.cpp
generator_LDADD=$(top_srcdir)/lib/libsimplecfd.la
generator_LDFLAGS=-lblas

if HAVE_PAPI
generator_SOURCES += event_set.cpp
generator_LDFLAGS += ${PAPI_LDFLAGS} ${PAPI_LDLIBS}
endif

