#include <simple_cfd/mesh.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/mesh_entity.hpp>
#include <string>

class Tester
{
private:
  const double epsilon;
  typedef cfd::TriangularCell cell_type;
  typedef cell_type::vertex_type vertex_type;

  void assertTrue(const bool b);
  void assertFalse(const bool b);
  void assertEqual(const double d, const double d2);

  void testQuadrature();  
  void testLinearBasis();
  void testQuadraticBasis();
  void testRisingFactorial();

  template<typename basis_t>
  void testBasis(const std::string& name)
  {
    std::cout << "Testing " << name << " basis values..." << std::endl;
  
    const double width = 20.0;
    const double height = 10.0;
    
    typedef typename basis_t::cell_type cell_type;
    cfd::TriangularMeshBuilder meshBuilder(width, height, 0.14);
    cfd::mesh<cell_type> m(meshBuilder.buildMesh());
    const std::size_t dimension = m.getDimension();
    basis_t basis(m);
  
    for(typename cfd::mesh<cell_type>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      const int dofs = basis.space_dimension();
      std::map<vertex_type, double> quadrature(m.getQuadrature(*cellIter));
  
      for(std::map<vertex_type, double>::const_iterator wIter(quadrature.begin()); wIter!=quadrature.end(); ++wIter)
      {
        double sum = 0.0;
        for(int i=0; i<dofs; ++i)
        {
          sum += basis.evaluate_tensor(cellIter->getIndex(), i, wIter->first).toScalar();
        }
        assertEqual(1.0, sum);
      }
    }
  }
public:
  Tester();
  void run();
};
