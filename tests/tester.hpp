#include <simple_cfd/mesh.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <string>

class Tester
{
private:
  const double epsilon;
  typedef cfd::cell<cfd::triangle> cell_type;
  typedef cell_type::vertex_type vertex_type;

  void assertTrue(const bool b);
  void assertFalse(const bool b);
  void assertZero(const double d);

  void testQuadrature();  
  void testLinearBasis();
  void testQuadraticBasis();

  template<typename basis_t>
  void testBasis(const std::string& name)
  {
    std::cout << "Testing " << name << " basis values..." << std::endl;
  
    const double width = 20.0;
    const double height = 10.0;
    
    typedef typename basis_t::cell_type cell_type;
    cfd::TriangularMeshBuilder meshBuilder(width, height, 0.14);
    cfd::mesh<cell_type> m(meshBuilder.buildMesh());
    basis_t basis(m);
    const std::map<cfd::cell_id, cell_type> cells(m.getCells());
  
    for(typename std::map<cfd::cell_id, cell_type>::const_iterator cellIter(cells.begin()); cellIter != cells.end(); ++cellIter)
    {
      const int dofs = basis.space_dimension();
      std::map<vertex_type, double> quadrature(cellIter->second.getQuadrature(m.getGeometry()));
  
      for(std::map<vertex_type, double>::const_iterator wIter(quadrature.begin()); wIter!=quadrature.end(); ++wIter)
      {
        double sum = 0.0;
        for(int i=0; i<dofs; ++i)
        {
          sum += basis.evaluate_tensor(cellIter->second, i, wIter->first).toScalar();
        }
        assertZero(1.0 - sum);
      }
    }
  }
public:
  Tester();
  void run();
};
