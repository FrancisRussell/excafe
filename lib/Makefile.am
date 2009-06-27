lib_LTLIBRARIES = libsimplecfd.la

CXX=${MPICXX}
libsimplecfd_la_SOURCES = mesh_connectivity.cpp mesh_topology.cpp mesh_entity_iterator_local.cpp petsc_manager.cpp \
triangular_mesh_builder.cpp numeric/matrix.cpp numeric/vector.cpp numeric/sparsity_pattern.cpp numeric/solver.cpp \
triangular_cell.cpp numeric/math_utilities.cpp numeric/polynomial.cpp
libsimplecfd_la_CXXFLAGS = -I/usr/include/petsc -I$(top_srcdir)/include
libsimplecfd_la_LDFLAGS = -ltriangle -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lX11 \
                          -lm -lblas -llapack 
