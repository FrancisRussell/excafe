#include <cstddef>
#include <sstream>
#include <boost/format.hpp>
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
#include <simple_cfd/exception.hpp>
#include <simple_cfd/boundary_condition_list.hpp>
#include <simple_cfd/boundary_condition_trivial.hpp>

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

  BoundaryCondition velocityConditions;

  SolveOperation coupledSolve;

  BoundaryCondition buildBoundaryConditions()
  {
    const Tensor<dimension> zero(1);

    Tensor<dimension> inflow(1);
    inflow(0) = 5.0;

    BoundaryConditionList<dimension> velocityConditionList(1);
    velocityConditionList.add(BoundaryConditionTrivial<dimension>(1, zero));
    velocityConditionList.add(BoundaryConditionTrivial<dimension>(3, zero));
    velocityConditionList.add(BoundaryConditionTrivial<dimension>(4, inflow));
    velocityConditionList.add(BoundaryConditionTrivial<dimension>(5, zero));

    return scenario.addBoundaryCondition(velocitySpace, velocityConditionList);
  }

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

    BoundaryConditionList<dimension> velocityConditionList(1);
    velocityConditions = buildBoundaryConditions();

    coupledSolve = constructCoupledSolver();
  }

  SolveOperation constructCoupledSolver()
  {
    using namespace forms;

    SolveOperation s = scenario.newSolveOperation();

    Scalar theta = 0.5;
    Scalar k  = 0.01;
    Scalar kinematic_viscosity = 1.0/250;

    Operator nonLinearRhs(velocitySpace, velocitySpace);
    nonLinearRhs =
      B(velocity, velocity)*dx +
      B(-(1.0-theta) * k * kinematic_viscosity * grad(velocity), grad(velocity))*dx +
      B((1.0-theta) * k * kinematic_viscosity * inner(grad(velocity), n), velocity)*ds +
      B(-(1.0-theta)*k * inner(velocityField, grad(velocity)), velocity)*dx;

    Field velocityRhs = nonLinearRhs * velocityField;
    Field load(project(velocityRhs, coupledSpace));

    TemporalIndex i;
    IndexedField unknownGuess(i);

    unknownGuess[-1] = project(velocityField, coupledSpace) + project(pressureField, coupledSpace);
    const forms::BilinearFormIntegralSum lhsForm = 
      B(velocity, velocity)*dx +
      B(theta * k * kinematic_viscosity * grad(velocity), grad(velocity))*dx +
      B(theta * k * kinematic_viscosity * -1.0 * inner(grad(velocity), n), velocity)*ds +
      B(-1.0 * k * pressure, div(velocity))*dx +
      B(div(velocity), pressure)*dx +
      B(theta * k * inner(project(unknownGuess[i-1], velocitySpace), grad(velocity)), velocity)*dx;

    LinearSystem system = assembleGalerkinSystem(coupledSpace, lhsForm, load, velocityConditions, unknownGuess[i-1]);
    Operator linearisedSystem = system.getConstrainedSystem();
    unknownGuess[i] = system.getSolution();

    Scalar residual = ((linearisedSystem * unknownGuess[i-1]) - system.getConstrainedLoad()).two_norm();

    i.setTermination(residual < 1e-3);

    s.setNewValue(velocityField, project(unknownGuess[final-1], velocitySpace));
    s.setNewValue(pressureField, project(unknownGuess[final-1], pressureSpace));

    s.finish();
    return s;
  }

  void step()
  {
    coupledSolve.execute();
  }

  void outputFieldsToFile(const std::string& filename)
  {
    scenario.outputFieldsToFile(filename);
  }
};

int main(int argc, char** argv)
{
  try
  {
    PETScManager::instance().init(argc, argv);
  
    static const std::size_t dimension = TriangularMeshBuilder::cell_dimension;
  
    TriangularMeshBuilder meshBuilder(3.0, 1.0, 1.0/900.0);
    meshBuilder.addPolygon(Polygon(vertex<2>(0.5, 0.5), 16, 0.148, 0), 5);
    Mesh<dimension> mesh(meshBuilder.buildMesh());
  
    NavierStokesSolver<dimension> solver(mesh);
  
    for(int i=0; i<6000; ++i)
    {
      std::cout << "Starting timestep " << i << "..." << std::endl;
      solver.step();
      std::ostringstream filename;
      filename << "./navier_stokes_" << boost::format("%|04|") % i << ".vtk";
      solver.outputFieldsToFile(filename.str());
    }
  }
  catch(const CFDException& e)
  {
    std::cerr << "A simple_cfd specific exception was generated: " << std::endl;
    std::cerr << e.what() << std::endl;
  }
}
