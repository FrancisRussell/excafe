#include <cstddef>
#include <simple_cfd/capture/scenario.hpp>
#include <simple_cfd/capture/solve_operation.hpp>
#include <simple_cfd/petsc_manager.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/lagrange_triangle_linear.hpp>
#include <simple_cfd/lagrange_triangle_quadratic.hpp>
#include <simple_cfd/capture/scenario.hpp>
#include <simple_cfd/capture/fields/function_space.hpp>
#include <simple_cfd/capture/fields/field.hpp>
#include <simple_cfd/mesh.hpp>

using namespace cfd;

template<std::size_t D>
class NavierStokesSolver
{
private:
  static const std::size_t dimension = D;
  Mesh<dimension> mesh;
  Scenario<dimension> scenario;

  Element velocity;
  Element pressure;

  FunctionSpace velocitySpace;
  FunctionSpace pressureSpace;
  FunctionSpace coupledSpace;

  NamedField velocityField;
  NamedField pressureField;

  SolveOperation coupledSolve;

public:
  NavierStokesSolver(Mesh<dimension>& _mesh) : mesh(_mesh), scenario(mesh)
  {
    velocity = scenario.addElement(new LagrangeTriangleQuadratic<1>());
    pressure = scenario.addElement(new LagrangeTriangleLinear<0>());

    velocitySpace = scenario.defineFunctionSpace(velocity, mesh);
    pressureSpace = scenario.defineFunctionSpace(pressure, mesh);
    coupledSpace = velocitySpace + pressureSpace;

    velocityField = scenario.defineNamedField("velocity", velocitySpace);
    pressureField = scenario.defineNamedField("pressure", pressureSpace);

    coupledSolve = constructCoupledSolver();
  }

  SolveOperation constructCoupledSolver()
  {
    SolveOperation s = scenario.newSolveOperation();
    return s;
  }
};

int main(int argc, char** argv)
{
  PETScManager::instance().init(argc, argv);

  static const std::size_t dimension = TriangularMeshBuilder::cell_dimension;

  TriangularMeshBuilder meshBuilder(3.0, 1.0, 1.0/900.0);
  meshBuilder.addPolygon(Polygon(vertex<2>(0.5, 0.5), 16, 0.148, 0), 5);
  Mesh<dimension> mesh(meshBuilder.buildMesh());

  NavierStokesSolver<dimension> solver(mesh);
}
