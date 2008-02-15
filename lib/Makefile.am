lib_LTLIBRARIES = libsimplecfd.la

CXX=mpicxx
libsimplecfd_la_SOURCES = petsc_manager.cpp
libsimplecfd_la_CXXFLAGS = -I/usr/include/petsc -Wno-long-long
libsimplecfd_la_LDFLAGS = -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lX11 -lm -lblas -llapack 
