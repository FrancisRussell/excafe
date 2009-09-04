#include <cstddef>
#include <simple_cfd/capture/scenario.hpp>
#include <simple_cfd/capture/solve_operation.hpp>
#include <simple_cfd/petsc_manager.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/lagrange_triangle_linear.hpp>
#include <simple_cfd/lagrange_triangle_quadratic.hpp>
#include <simple_cfd/capture/scenario.hpp>
#include <simple_cfd/capture/fields/fields.hpp>
#include <simple_cfd/capture/forms/forms.hpp>
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
    using namespace forms;

    SolveOperation s = scenario.newSolveOperation();

    Scalar theta = 0.5;
    Scalar k  = 0.01;
    Scalar kinematic_viscosity = 1.0/250;

    Operator systemMatrix(coupledSpace, coupledSpace);
    systemMatrix = 
      B(velocity, velocity)*dx +
      B(theta * k * kinematic_viscosity * grad(velocity), grad(velocity))*dx +
      B(theta * k * kinematic_viscosity * -1.0 * inner(grad(velocity), n), velocity)*ds +
      B(-1.0 * k * pressure, div(velocity))*dx +
      B(div(velocity), pressure)*dx;

    Operator nonLinearRhs(velocitySpace, velocitySpace);
    nonLinearRhs =
      B(velocity, velocity)*dx +
      B(-(1.0-theta) * k * kinematic_viscosity * grad(velocity), grad(velocity))*dx +
      B((1.0-theta) * k * kinematic_viscosity * inner(grad(velocity), n), velocity)*ds +
      B(-(1.0-theta)*k * inner(velocityField, grad(velocity)), velocity)*dx;

    Field velocityRhs = nonLinearRhs * velocityField;
    Field load(coupledSpace);

    TemporalIndex n;
    IndexedField unknownGuess(n);

    unknownGuess[-1] = project(velocityField, coupledSpace) + project(pressureField, coupledSpace);
    Operator linearisedSystem = 
      systemMatrix + B(theta * k * inner(project(unknownGuess[n-1], velocitySpace), grad(velocity)), velocity)*dx;

    Scalar residual = ((linearisedSystem * unknownGuess[n-1]) - load).two_norm();

    // TODO: Boundary condition magic
    unknownGuess[n] = linear_solve(linearisedSystem, load);

    n.setTermination(residual < 1e-3);

    s.setNewValue(velocityField, project(unknownGuess[final-1], velocitySpace));
    s.setNewValue(pressureField, project(unknownGuess[final-1], pressureSpace));

    s.finish();
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
