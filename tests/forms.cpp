#include <excafe/capture/forms/forms.hpp>
#include <excafe/lagrange_triangle_linear.hpp>
#include <excafe/mesh.hpp>
#include <excafe/triangular_mesh_builder.hpp>
#include <excafe/dof_map.hpp>
#include <excafe/dof_map_builder.hpp>
#include <excafe/discrete_field.hpp>
#include <excafe/petsc_manager.hpp>
#include <excafe/numeric/tensor.hpp>

int main(int argc, char** argv)
{
  using namespace excafe;
  using namespace excafe::forms;

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
