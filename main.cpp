#include <mesh.hpp>
#include <dof_map_builder.hpp>
#include <lagrange_triangle_linear.hpp>
#include <lagrange_triangle_quadratic.hpp>

using namespace cfd;

int main(int argc, char* argv[])
{
  typedef cell<triangle> cell_type;
  mesh<cell_type> m(7,7);
  m.print();

  lagrange_triangle_linear pressure(m);
  lagrange_triangle_quadratic velocity_x(m);
  lagrange_triangle_quadratic velocity_y(m);

  dof_map_builder<cell_type> mapBuilder(m);
  mapBuilder.addFiniteElement(pressure);
  mapBuilder.addFiniteElement(velocity_x);
  mapBuilder.addFiniteElement(velocity_y);
  mapBuilder.handleCells(m.getCells());

  dof_map<cell_type> dofMap(mapBuilder.getDofMap());
  std::cout << "Size of dof map: " << dofMap.getMappingSize() << std::endl;
  std::cout << "Degrees of freedom: " << dofMap.getDegreesOfFreedom() << std::endl;
}
