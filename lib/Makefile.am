lib_LTLIBRARIES = libsimplecfd.la

CXX=${MPICXX}
libsimplecfd_la_SOURCES = mesh_connectivity.cpp mesh_topology.cpp mesh_entity_iterator_local.cpp petsc_manager.cpp \
triangular_mesh_builder.cpp numeric/matrix.cpp numeric/vector.cpp numeric/sparsity_pattern.cpp numeric/solver.cpp \
triangular_cell.cpp numeric/math_utilities.cpp \
numeric/quadrature.cpp capture/fields/operator.cpp capture/fields/scalar.cpp \
capture/fields/field.cpp capture/fields/function_space.cpp \
capture/fields/discrete_expr_container_builder_form_visitor.cpp capture/fields/discrete_expr.cpp \
capture/fields/operator_assembly.cpp capture/evaluation/discrete_expr_scoping.cpp \
capture/fields/temporal_index_set.cpp capture/fields/temporal_index_value.cpp \
capture/evaluation/evaluation_strategy.cpp capture/solve_operation.cpp capture/array/array_index.cpp \
capture/array/tensor_index.cpp

libsimplecfd_la_CXXFLAGS = -I/usr/include/petsc -I$(top_srcdir)/include
libsimplecfd_la_LDFLAGS = -ltriangle -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lX11 \
                          -lm -lblas -llapack 
