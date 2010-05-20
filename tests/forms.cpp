#include <simple_cfd/capture/forms/forms.hpp>
#include <simple_cfd/lagrange_triangle_linear.hpp>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/dof_map.hpp>
#include <simple_cfd/dof_map_builder.hpp>
#include <simple_cfd/discrete_field.hpp>
#include <simple_cfd/petsc_manager.hpp>
#include <simple_cfd/numeric/tensor.hpp>

int main(int argc, char** argv)
{
  using namespace cfd;
  using namespace cfd::forms;

  PETScManager::instance().init(argc, argv);

  // Build mesh
  typedef TriangularCell cell_type;
  static const std::size_t dimension = cell_type::dimension;

  TriangularMeshBuilder meshBuilder(1.0, 1.0, 0.25);
  Mesh<dimension> m(meshBuilder.buildMesh());

  // Define bases and dof maps
  LagrangeTriangleLinear<1> vectorBasis;
  DofMapBuilder<dimension> mapBuilder(m);
  mapBuilder.addFiniteElement(vectorBasis);
  DofMap<dimension> vectorFieldDofs = mapBuilder.getDofMap();

  DiscreteField<dimension> vectorField(vectorFieldDofs);

  grad(vectorField);
  div(vectorBasis) + div(vectorField);

  Tensor<3> rankTwoTensor(2);
  div(rankTwoTensor);

  B(div(rankTwoTensor), div(rankTwoTensor))*dx + B(div(rankTwoTensor), div(rankTwoTensor))*ds;
}
