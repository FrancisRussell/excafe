lib_LTLIBRARIES = libsimplecfd.la
noinst_LTLIBRARIES = libcse.la libsymbolic.la

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
util/type_info.cpp util/timer.cpp util/aligned_alloc.cpp numeric/factoriser.cpp \
numeric/factoriser_table.cpp codegen/ufc_integral_generator.cpp codegen/dynamic_cxx.cpp \
codegen/product.cpp util/apr_manager.cpp util/apr_pool.cpp \
capture/assembly/generic_symbol.cpp \
mp/integer.cpp mp/rational.cpp mp/float.cpp mp/cln_conversions.cpp

libcse_la_CXXFLAGS = -I$(top_srcdir)/include
libcse_la_SOURCES = cse/cube.cpp cse/sop.cpp cse/polynomial_index.cpp cse/sop_map.cpp cse/kcm.cpp

libsymbolic_la_CXXFLAGS = -I$(top_srcdir)/include
libsymbolic_la_SOURCES= symbolic/expr.cpp symbolic/symbol.cpp symbolic/float.cpp \
symbolic/basic.cpp symbolic/sum.cpp symbolic/product.cpp \
symbolic/expand_visitor.cpp symbolic/rational.cpp symbolic/group.cpp symbolic/collect_visitor.cpp \
symbolic/polynomial.cpp symbolic/abstract_basic.cpp \
symbolic/derivative_cache.cpp \
symbolic/abs.cpp symbolic/extracted_expressions.cpp

libsimplecfd_la_CXXFLAGS = -I/usr/include/petsc -I$(top_srcdir)/include
libsimplecfd_la_LIBADD = libcse.la libsymbolic.la
libsimplecfd_la_LDFLAGS = -ltriangle -lm ${PETSC_LIBS} ${GMP_LIBS} ${CLN_LIBS} ${GINAC_LIBS} ${APR1_LIBS}
