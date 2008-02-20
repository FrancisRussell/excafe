lib_LTLIBRARIES = libsimplecfd.la

CXX=mpicxx
libsimplecfd_la_SOURCES = petsc_manager.cpp numeric/matrix.cpp numeric/vector.cpp numeric/sparsity_pattern.cpp numeric/solver.cpp
libsimplecfd_la_CXXFLAGS = -I/usr/include/petsc -I$(top_srcdir)/include -Wno-long-long
libsimplecfd_la_LDFLAGS = -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lX11 -lm -lblas -llapack 
