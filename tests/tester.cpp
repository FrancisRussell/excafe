#include <tester.hpp>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/lagrange_triangle_linear.hpp>
#include <simple_cfd/lagrange_triangle_quadratic.hpp>
#include <map>
#include <iostream>
#include <cmath>

using namespace cfd;

Tester::Tester() : epsilon(1e-10)
{
}

void Tester::run() 
{
  testQuadrature();
  testLinearBasis();
  testQuadraticBasis();
}

void Tester::assertTrue(const bool b)
{
  if (!b)
  {
    std::cout << "assertTrue failed" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Tester::assertFalse(const bool b)
{
  if (b)
  {
    std::cout << "assertFalse failed" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Tester::assertZero(const double d)
{
  if (std::abs(d) >= epsilon)
  {
    std::cout << "assertZero failed" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void Tester::testQuadrature()
{
  std::cout << "Testing quadrature..." << std::endl;

  const double area = 50.0;
  const double width = 10.0;
  TriangularMeshBuilder meshBuilder(width, width, area);
  mesh<cell_type> m(meshBuilder.buildMesh());
  const std::size_t dimension = m.getDimension();

  for(mesh<cell_type>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
  {
    // Check the area is correct first
    assertZero(m.getArea(cellIter->getIndex()) - area);
    std::map<vertex_type, double> quadrature(m.getQuadrature(*cellIter));

    double accum = 0;
    for(std::map<vertex_type, double>::const_iterator wIter(quadrature.begin()); wIter!=quadrature.end(); ++wIter)
    {
      accum += wIter->second * m.getReferenceCell().getJacobian(m, *cellIter, wIter->first);
    }

    assertZero(accum - area);
  }
}

void Tester::testLinearBasis()
{
  testBasis< lagrange_triangle_linear<0> >("linear");
}

void Tester::testQuadraticBasis()
{
  testBasis< lagrange_triangle_quadratic<0> >("quadratic");
}
