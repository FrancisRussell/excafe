#include <simple_cfd/mesh.hpp>
#include <simple_cfd/stokes_system.hpp>
#include <simple_cfd/petsc_manager.hpp>
#include <iostream>
#include <sstream>

using namespace cfd;

int main(int argc, char* argv[])
{
  PETScManager::instance().init(argc, argv);
  typedef cell<triangle> cell_type;
  mesh<cell_type> m(3.0, 1.0, 90, 30);

  std::cout << "Constructing system..." << std::endl;
  stokes_system<cell_type> system(m);
  std::cout << "Assembling system..." << std::endl;
  //std::cout << "Applying boundary conditions..." << std::endl;
  //system.applyBoundaryConditions();

  //std::cout << "Starting solver..." << std::endl;
  //system.solve();
  for(int i=0; i<120; ++i)
  {
    std::cout << "Starting timestep: " << i << std::endl;
    system.timeDependentAssembleAndSolve();
    std::stringstream filename;
    filename << "./steady_stokes_" << i << ".vtk";
    system.outputToFile(filename.str());
  }
}
