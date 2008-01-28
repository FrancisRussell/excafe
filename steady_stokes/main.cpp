#include <simple_cfd/mesh.hpp>
#include <simple_cfd/stokes_system.hpp>
#include <iostream>

using namespace cfd;

int main(int argc, char* argv[])
{
  typedef cell<triangle> cell_type;
  mesh<cell_type> m(50,50);

  std::cout << "Constructing system..." << std::endl;
  stokes_system<cell_type> system(m);
  std::cout << "Assembling system..." << std::endl;
  system.assemble();
  std::cout << "Applying boundary conditions..." << std::endl;
  system.applyBoundaryConditions();

  std::cout << "Starting solver..." << std::endl;
  system.solve();
  system.outputToFile("./steady_stokes.vtk");
}
