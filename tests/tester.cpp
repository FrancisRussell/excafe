#include <tester.hpp>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/lagrange_triangle_linear.hpp>
#include <simple_cfd/lagrange_triangle_quadratic.hpp>
#include <simple_cfd/numeric/math_utilities.hpp>
#include <simple_cfd/numeric/quadrature.hpp>
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
  testTriangleQuadrature();
  testLinearBasis();
  testQuadraticBasis();
  testLinearBasisDofs();
  testQuadraticBasisDofs();
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

void Tester::assertEqual(const double d, const double d2)
{
  if (std::abs(d-d2) >= epsilon)
  {
    std::cout << "assertEqual failed: " << d << "!=" << d2 << std::endl;
    exit(EXIT_FAILURE);
  }
}


void Tester::testTriangleQuadrature()
{
  std::cout << "Testing triangle quadrature..." << std::endl;

  const double area = 50.0;
  const double width = 10.0;
  TriangularMeshBuilder meshBuilder(width, width, area);
  mesh<cell_type::dimension> m(meshBuilder.buildMesh());
  const std::size_t dimension = m.getDimension();

  cfd::QuadraturePoints<2> quadrature(m.getReferenceCell().getQuadrature(5));
  const cfd::MeshEntity localCell(dimension, 0);

  for(mesh<cell_type::dimension>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
  {
    // Check the area is correct first
    assertEqual(m.getArea(cellIter->getIndex()), area);
    const double jacobian = m.getJacobian(cellIter->getIndex(), vertex_type(0.0, 0.0));
    double accum = 0;

    for(cfd::QuadraturePoints<2>::const_iterator wIter(quadrature.begin(localCell)); wIter!=quadrature.end(localCell); ++wIter)
      accum += wIter->second * jacobian;

    assertEqual(accum, area);
  }
}

void Tester::testQuadrature()
{

  Quadrature quadrature;

  std::cout << "Testing Gauss-Legendre quadrature..." << std::endl;
  testQuadrature(quadrature.getGauss(5));

  std::cout << "Testing Gauss-Radau-Legendre quadrature..." << std::endl;
  testQuadrature(quadrature.getGaussRadau(5));

  std::cout << "Testing Gauss-Lobatto-Legendre quadrature..." << std::endl;
  testQuadrature(quadrature.getGaussLobatto(5));
}

void Tester::testQuadrature(const std::map<double, double>& q)
{
  double sum = 0.0;

  for(std::map<double, double>::const_iterator qIter(q.begin()); qIter!=q.end(); ++qIter)
    sum += qIter->second;

  assertEqual(sum, 2.0);
}

void Tester::testLinearBasis()
{
  testBasis< lagrange_triangle_linear<0> >("linear");
}

void Tester::testQuadraticBasis()
{
  testBasis< lagrange_triangle_quadratic<0> >("quadratic");
}

void Tester::testLinearBasisDofs()
{
  testBasisDofs< lagrange_triangle_linear<0> >("linear");
}

void Tester::testQuadraticBasisDofs()
{
  testBasisDofs< lagrange_triangle_quadratic<0> >("quadratic");
}
