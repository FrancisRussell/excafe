bin_PROGRAMS = main
main_SOURCES = main.cpp 
main_CXXFLAGS = -I.
main_LDFLAGS = -ldolfin

STEADY_STOKES_GENERATED = SteadyStokes.h
BUILT_SOURCES = ${STEADY_STOKES_GENERATED}
CLEANFILES = ${STEADY_STOKES_GENERATED}
${STEADY_STOKES_GENERATED}: SteadyStokes.ufl
	ffc -l dolfin -r quadrature SteadyStokes.ufl
