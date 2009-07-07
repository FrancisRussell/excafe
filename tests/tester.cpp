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
  testRisingFactorial();
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
  mesh<cell_type> m(meshBuilder.buildMesh());
  const std::size_t dimension = m.getDimension();

  for(mesh<cell_type>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
  {
    // Check the area is correct first
    assertEqual(m.getArea(cellIter->getIndex()), area);
    std::map<vertex_type, double> quadrature(m.getQuadrature(*cellIter));

    double accum = 0;
    for(std::map<vertex_type, double>::const_iterator wIter(quadrature.begin()); wIter!=quadrature.end(); ++wIter)
    {
      accum += wIter->second;
    }

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

void Tester::testRisingFactorial()
{
  std::cout << "Testing rising factorial..." << std::endl;

  double x = 5.0;
  assertEqual(MathUtilities::rising_factorial(x, 0), 1.0);
  assertEqual(MathUtilities::rising_factorial(x, 1), x);
  assertEqual(MathUtilities::rising_factorial(x, 2), x*x + x);
  assertEqual(MathUtilities::rising_factorial(x, 3), x*x*x + 3.0*x*x + 2.0*x);
  assertEqual(MathUtilities::rising_factorial(x, 4), x*x*x*x + 6.0*x*x*x + 11.0*x*x + 6.0*x);
}
