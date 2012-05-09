#include <cstddef>
#include <sstream>
#include <boost/format.hpp>
#include <excafe/capture/scenario.hpp>
#include <excafe/capture/solve_operation.hpp>
#include <excafe/petsc_manager.hpp>
#include <excafe/triangular_mesh_builder.hpp>
#include <excafe/lagrange_triangle.hpp>
#include <excafe/capture/scenario.hpp>
#include <excafe/capture/fields/fields.hpp>
#include <excafe/capture/forms/forms.hpp>
#include <excafe/mesh.hpp>
#include <excafe/exception.hpp>
#include <excafe/boundary_condition_list.hpp>
#include <excafe/boundary_condition_trivial.hpp>

using namespace excafe;

template<std::size_t D>
class HeatSolver
{
private:
  static const std::size_t dimension = D;
  Mesh<dimension> mesh;
  Scenario<dimension> scenario;

  Element temperature;
  FunctionSpace temperatureSpace;
  NamedField temperatureField;

  BoundaryCondition boundaryConditions;

  SolveOperation solve;

  BoundaryCondition buildBoundaryConditions()
  {
    Tensor<dimension> zero(0);
    zero  = 0;

    Tensor<dimension> source(0);
    source = 1.0;

    BoundaryConditionList<dimension> boundaryConditionList(0);
    boundaryConditionList.add(BoundaryConditionTrivial<dimension>(1, source));
    boundaryConditionList.add(BoundaryConditionTrivial<dimension>(2, zero));
    boundaryConditionList.add(BoundaryConditionTrivial<dimension>(3, zero));
    boundaryConditionList.add(BoundaryConditionTrivial<dimension>(4, zero));
    boundaryConditionList.add(BoundaryConditionTrivial<dimension>(5, source));

    return scenario.addBoundaryCondition(temperatureSpace, boundaryConditionList);
  }

public:
  HeatSolver(Mesh<dimension>& _mesh) : mesh(_mesh), scenario(mesh)
  {
    temperature = scenario.addElement(new LagrangeTriangle<0>(1));
    temperatureSpace = scenario.defineFunctionSpace(temperature, mesh);
    temperatureField = scenario.defineNamedField("temperature", temperatureSpace);

    BoundaryConditionList<dimension> boundaryConditionList(0);
    boundaryConditions = buildBoundaryConditions();

    solve = constructSolve();
  }

  SolveOperation constructSolve()
  {
    using namespace forms;

    SolveOperation s = scenario.newSolveOperation();

    Scalar c = 1e-4;
    Scalar k  = 10.0;

    Operator massMatrix(temperatureSpace, temperatureSpace);
    massMatrix = B(temperature, temperature)*dx;

    const forms::BilinearFormIntegralSum lhsForm = 
      B(temperature, temperature)*dx + B(k*c*grad(temperature), grad(temperature))*dx;

    Field rhs = massMatrix * temperatureField;

    LinearSystem system = assembleGalerkinSystem(temperatureSpace, lhsForm, rhs, boundaryConditions, temperatureField);
    s.setNewValue(temperatureField, system.getSolution());

    s.finish();
    return s;
  }

  void step()
  {
    solve.execute();
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
    static const double maxCellArea = 1.0/5000;
    static const double polySize = 0.1;
    static const std::size_t polyEdges = 4;
    static const double polyLabel = 5;
    static const double polyRotation = 0.0;

    TriangularMeshBuilder meshBuilder(1.0, 1.0, maxCellArea);
    const Polygon poly(vertex<2>(0.5, 0.5), polyEdges, polySize, polyRotation);
    meshBuilder.addPolygon(poly, polyLabel);
    Mesh<dimension> mesh(meshBuilder.buildMesh());
  
    HeatSolver<dimension> solver(mesh);
  
    for(int i=0; i<200; ++i)
    {
      std::cout << "Starting timestep " << i << "..." << std::endl;
      solver.step();
      std::ostringstream filename;
      filename << "./heat_" << boost::format("%|04|") % i << ".vtk";
      solver.outputFieldsToFile(filename.str());
    }
  }
  catch(const CFDException& e)
  {
    std::cerr << "A excafe specific exception was generated: " << std::endl;
    std::cerr << e.what() << std::endl;
  }
}
