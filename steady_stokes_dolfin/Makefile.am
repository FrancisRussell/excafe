bin_PROGRAMS = main

main_SOURCES = main.cpp 
main_CXXFLAGS = ${LIBXML2_CFLAGS}
main_LDFLAGS = -ldolfin ${LIBXML2_LIBS}

STEADY_STOKES_GENERATED = SteadyStokes.h
BUILT_SOURCES = ${STEADY_STOKES_GENERATED}
CLEANFILES = ${STEADY_STOKES_GENERATED}
${STEADY_STOKES_GENERATED}: SteadyStokes.ufl
	ffc -O -l dolfin -r quadrature SteadyStokes.ufl
