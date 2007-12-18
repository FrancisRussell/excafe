#include <simple_cfd/mesh.hpp>
#include <simple_cfd/stokes_system.hpp>

using namespace cfd;

int main(int argc, char* argv[])
{
  typedef cell<triangle> cell_type;
  mesh<cell_type> m(7,7);
  m.print();

  stokes_system<cell_type> system(m);
  system.assemble();
  system.applyBoundaryConditions();

  system.solve();
}
