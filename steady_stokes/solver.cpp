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

class BoundaryConditionHack
{
private:
  static const std::size_t dimension = 2;
  typedef DofMap<dimension>::dof_t dof_t;
  typedef vertex<dimension> vertex_type;

  Mesh<dimension>& m;
  FiniteElement<dimension>& velocity;
  FiniteElement<dimension>& pressure;

  enum Location
  {
    TOP_EDGE,
    BOTTOM_EDGE,
    LEFT_EDGE,
    RIGHT_EDGE,
    BODY
  };

  Location getLocation(const vertex_type& v) const
  {
    // FIXME: make me into a function that doesn't depend on the specifics of the mesh generation

    Location location = BODY;

    if (v[0] == 3.0)
      location = RIGHT_EDGE;

    if (v[1] == 0.0)
      location = BOTTOM_EDGE;

    if (v[1] == 1.0)
      location = TOP_EDGE;

    if (v[0] == 0.0)
      location = LEFT_EDGE;

    return location;
  }

  void applyEdgeVelocityBoundaryConditions(DiscreteOperator<dimension>* const stiffness_matrix, 
    DiscreteField<dimension>* const load_vector) const
  {
    const unsigned velocitySpaceDimension = velocity.spaceDimension();

    for(Mesh<dimension>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      for(unsigned dof=0; dof<velocitySpaceDimension; ++dof)
      {
        const dof_t velocity_globalDof = dof_t(&velocity, cellIter->getIndex(), dof);
        const vertex_type dofLocation = velocity.getDofCoordinateGlobal(m, cellIter->getIndex(), dof);
        const bool isXDof = velocity.getTensorIndex(m, cellIter->getIndex(), dof) == 0;
        const Location location = getLocation(dofLocation);

        if (location == LEFT_EDGE)
        {
          if (stiffness_matrix != NULL)
            stiffness_matrix->zeroRow(velocity_globalDof, 1.0);

          // Set x velocity to same value on inflow and outflow boundary, and y velocity to 0
          const double rhs = isXDof ? 5.0 : 0.0;

          if (load_vector != NULL)
            load_vector->setValues(1, &velocity_globalDof, &rhs);

        }
        else if ((location == TOP_EDGE || location == BOTTOM_EDGE) && !isXDof)
        {
          if (stiffness_matrix != NULL)
            stiffness_matrix->zeroRow(velocity_globalDof, 1.0);

          const double rhs = 0.0;

          if (load_vector != NULL)
            load_vector->setValues(1, &velocity_globalDof, &rhs);
        }
      }
    }

    if (stiffness_matrix != NULL)
      stiffness_matrix->assemble();

    if (load_vector != NULL)
      load_vector->assemble();
  }

  void applyCylinderVelocityBoundaryConditions(DiscreteOperator<dimension>* const stiffness_matrix, 
    DiscreteField<dimension>* const load_vector) const
  {
    const unsigned velocitySpaceDimension = velocity.spaceDimension();
    const vertex_type centre(0.5, 0.5);
    const double radius = 0.15;

    for(Mesh<dimension>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      for(unsigned dof=0; dof<velocitySpaceDimension; ++dof)
      {
        const vertex_type dofLocation = velocity.getDofCoordinateGlobal(m, cellIter->getIndex(), dof);
        const vertex_type offset = dofLocation - centre;

        if((offset[0] * offset[0] + offset[1] * offset[1]) < radius * radius)
        {
          const dof_t velocity_globalDof = dof_t(&velocity, cellIter->getIndex(), dof);

          if (stiffness_matrix != NULL)
            stiffness_matrix->zeroRow(velocity_globalDof, 1.0);

          const double rhs = 0.0;

          if (load_vector != NULL)
            load_vector->setValues(1, &velocity_globalDof, &rhs);
        }
      }
    }

    if (stiffness_matrix != NULL)
      stiffness_matrix->assemble();

    if (load_vector != NULL)
      load_vector->assemble();
  }

public:
  BoundaryConditionHack(Mesh<dimension>& _m, FiniteElement<dimension>& _velocity, FiniteElement<dimension>& _pressure) :
   m(_m), velocity(_velocity), pressure(_pressure)
  {
    assert(pressure.getRank() == 0);
    assert(velocity.getRank() == 1);
  }

  void operator()(DiscreteOperator<dimension>& stiffness_matrix) const
  {
    applyEdgeVelocityBoundaryConditions(&stiffness_matrix, NULL);
    applyCylinderVelocityBoundaryConditions(&stiffness_matrix, NULL);
  }

  void operator()(DiscreteField<dimension>& load_vector) const
  {
    applyEdgeVelocityBoundaryConditions(NULL, &load_vector);
    applyCylinderVelocityBoundaryConditions(NULL, &load_vector);
  }
};

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

    const BoundaryConditionHack bc(mesh, scenario.getElement(velocity), scenario.getElement(pressure));

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

    LinearSystem system = assembleGalerkinSystem(coupledSpace, lhsForm, load, velocityConditions);
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
      solver.step();
      std::stringstream filename;
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
