#include <cstddef>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <simple_cfd/capture/scenario.hpp>
#include <simple_cfd/capture/solve_operation.hpp>
#include <simple_cfd/petsc_manager.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/lagrange_triangle.hpp>
#include <simple_cfd/capture/scenario.hpp>
#include <simple_cfd/capture/fields/fields.hpp>
#include <simple_cfd/capture/forms/forms.hpp>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/boundary_condition_list.hpp>
#include <simple_cfd/boundary_condition_trivial.hpp>

using namespace cfd;

template<std::size_t D>
class MassMatrixGenerator
{
private:
  static const std::size_t dimension = D;
  Mesh<dimension> mesh;
  Scenario<dimension> scenario;
  int nf;
  int p;
  int q;

  Element elementBasis;
  FunctionSpace elementSpace;

  Element coefficientBasis;
  FunctionSpace coefficientSpace;

  NamedField f;
  NamedField g;
  NamedField h;
  NamedField i;

public:
  MassMatrixGenerator(Mesh<dimension>& _mesh, const int _nf, const int _p, const int _q) :
    mesh(_mesh), scenario(mesh), nf(_nf), p(_p), q(_q)
  {
    elementBasis = scenario.addElement(new LagrangeTriangle<0>(q));
    elementSpace = scenario.defineFunctionSpace(elementBasis, mesh);

    coefficientBasis = scenario.addElement(new LagrangeTriangle<0>(p));
    coefficientSpace = scenario.defineFunctionSpace(coefficientBasis, mesh);

    f = scenario.defineNamedField("f", coefficientSpace);
    g = scenario.defineNamedField("g", coefficientSpace);
    h = scenario.defineNamedField("h", coefficientSpace);
    i = scenario.defineNamedField("i", coefficientSpace);
  }

  forms::BilinearFormIntegralSum constructForm()
  {
    using namespace forms;

    LinearForm trial = elementBasis;
    if (nf > 0) trial = trial*f;
    if (nf > 1) trial = trial*g;
    if (nf > 2) trial = trial*h;
    if (nf > 3) trial = trial*i;
    if (nf > 4) CFD_EXCEPTION("Cannot generate code for more than 4 coefficient functions.");
    
    const BilinearFormIntegralSum form = B(trial, elementBasis)*dx;
    return form;
  }

  void dumpUFL(std::ostream& out)
  {
    const forms::BilinearFormIntegralSum form = constructForm();
    scenario.writeUFCCellIntegral(out, form);
  }
};

int main(int argc, char** argv)
{
  if (argc != 5)
  {
    std::cerr << "Usage: generator num_functions coefficient_degree element_degree outfile" << std::endl;
    exit(1);
  }
  else
  {
    const int nf = boost::lexical_cast<int>(argv[1]);
    const int p = boost::lexical_cast<int>(argv[2]);
    const int q = boost::lexical_cast<int>(argv[3]);
    const std::string filename = boost::lexical_cast<std::string>(argv[4]);

    try
    {
      PETScManager::instance().init(argc, argv);
    
      static const std::size_t dimension = TriangularMeshBuilder::cell_dimension;
      static const double maxCellArea = 1.0/5000;
  
      TriangularMeshBuilder meshBuilder(1.0, 1.0, maxCellArea);
      Mesh<dimension> mesh(meshBuilder.buildMesh());
    
      MassMatrixGenerator<dimension> generator(mesh, nf, p, q);
      std::ofstream out(filename.c_str());
      generator.dumpUFL(out);
      out.close();
    }
    catch(const CFDException& e)
    {
      std::cerr << "A simple_cfd specific exception was generated: " << std::endl;
      std::cerr << e.what() << std::endl;
    }
  }
}
