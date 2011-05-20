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
capture/evaluation/evaluation_strategy.cpp capture/solve_operation.cpp \
capture/assembly/scalar_placeholder.cpp capture/assembly/position_component.cpp \
capture/assembly/cell_vertex_component.cpp capture/assembly/basis_coefficient.cpp \
capture/assembly/scalar_access.cpp numeric/tensor_size.cpp capture/indices/propagation_rules.cpp \
cse/cube.cpp cse/sop.cpp cse/polynomial_index.cpp \
symbolic/expr.cpp symbolic/symbol.cpp symbolic/float.cpp \
util/type_info.cpp symbolic/basic.cpp symbolic/sum.cpp symbolic/product.cpp \
symbolic/expand_visitor.cpp symbolic/rational.cpp symbolic/group.cpp symbolic/collect_visitor.cpp \
symbolic/polynomial.cpp symbolic/abstract_basic.cpp cse/sop_map.cpp numeric/factoriser.cpp \
numeric/factoriser_table.cpp codegen/ufc_integral_generator.cpp codegen/dynamic_cxx.cpp \
util/apr_manager.cpp util/apr_pool.cpp symbolic/derivative_cache.cpp

libsimplecfd_la_CXXFLAGS = -I/usr/include/petsc -I$(top_srcdir)/include
libsimplecfd_la_LDFLAGS = -ltriangle -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc -lX11 \
                          -lm -lblas -llapack ${CLN_LIBS} ${GINAC_LIBS} ${APR1_LIBS}
