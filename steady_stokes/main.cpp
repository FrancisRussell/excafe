#include <simple_cfd/mesh.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/stokes_system.hpp>
#include <simple_cfd/petsc_manager.hpp>
#include <simple_cfd/subdomain.hpp>
#include <simple_cfd/function.hpp>
#include <simple_cfd/numeric/tensor.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <sstream>

using namespace cfd;

int main(int argc, char* argv[])
{
  PETScManager::instance().init(argc, argv);
  typedef TriangularCell cell_type;
  TriangularMeshBuilder meshBuilder(3.0, 1.0, 1.0/900.0);
  meshBuilder.addPolygon(Polygon(vertex<2>(1.0, 0.5), 16, 0.148, 0), 5);
  mesh<cell_type::dimension> m(meshBuilder.buildMesh());

  std::cout << "Constructing system..." << std::endl;
  stokes_system<cell_type> system(m);

  for(int i=0; i<6000; ++i)
  {
    std::cout << "Starting timestep: " << i << std::endl;
    system.coupledSolve();
    std::stringstream filename;
    filename << "./steady_stokes_" << boost::format("%|04|") % i << ".vtk";
    system.outputToFile(filename.str());
  }
}
