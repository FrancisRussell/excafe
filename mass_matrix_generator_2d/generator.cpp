#include <cstddef>
#include <sstream>
#include <boost/lexical_cast.hpp>
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
class MassMatrixGenerator
{
private:
  static const std::size_t dimension = D;
  Mesh<dimension> mesh;
  Scenario<dimension> scenario;
  bool laplacian;
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
  MassMatrixGenerator(Mesh<dimension>& _mesh, const bool _laplacian, const int _nf, const int _p, const int _q) :
    mesh(_mesh), scenario(mesh), laplacian(_laplacian), nf(_nf), p(_p), q(_q)
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

    const LinearForm transformedBasis = (laplacian ? grad(elementBasis) : elementBasis);
    LinearForm trial = transformedBasis;
    if (nf > 0) trial = trial*f;
    if (nf > 1) trial = trial*g;
    if (nf > 2) trial = trial*h;
    if (nf > 3) trial = trial*i;
    if (nf > 4) CFD_EXCEPTION("Cannot generate code for more than 4 coefficient functions.");
    
    const BilinearFormIntegralSum form = B(trial, transformedBasis)*dx;
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
  if (argc != 6)
  {
    std::cerr << "Usage: generator (mass_matrix|laplacian) num_functions coefficient_degree element_degree outfile" << std::endl;
    exit(1);
  }
  else
  {
    const std::string matType(argv[1]);
    bool laplacian;

    if (matType == "mass_matrix")
      laplacian = false;
    else if (matType == "laplacian")
      laplacian = true;
    else
      CFD_EXCEPTION("First argument must be mass_matrix or laplacian.");

    const int nf = boost::lexical_cast<int>(argv[2]);
    const int p = boost::lexical_cast<int>(argv[3]);
    const int q = boost::lexical_cast<int>(argv[4]);
    const std::string filename(argv[5]);

    try
    {
      PETScManager::instance().init(argc, argv);
    
      static const std::size_t dimension = TriangularMeshBuilder::cell_dimension;
      static const double maxCellArea = 1.0/5000;
  
      TriangularMeshBuilder meshBuilder(1.0, 1.0, maxCellArea);
      Mesh<dimension> mesh(meshBuilder.buildMesh());
    
      MassMatrixGenerator<dimension> generator(mesh, laplacian, nf, p, q);
      std::ofstream out(filename.c_str());
      generator.dumpUFL(out);
      out.close();
    }
    catch(const CFDException& e)
    {
      std::cerr << "A excafe specific exception was generated: " << std::endl;
      std::cerr << e.what() << std::endl;
    }
  }
}
