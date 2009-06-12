#include <simple_cfd/mesh.hpp>
#include <simple_cfd/triangular_mesh_builder.hpp>
#include <simple_cfd/dof_map_builder.hpp>
#include <simple_cfd/mesh_entity.hpp>
#include <string>
#include <vector>

class Tester
{
private:
  const double epsilon;
  typedef cfd::TriangularCell cell_type;
  typedef cell_type::vertex_type vertex_type;

  void assertTrue(const bool b);
  void assertFalse(const bool b);
  void assertZero(const double d);

  void testQuadrature();  
  void testLinearBasis();
  void testQuadraticBasis();
  void testLinearBasisDofs();
  void testQuadraticBasisDofs();

  template<typename basis_t>
  void testBasis(const std::string& name)
  {
    std::cout << "Testing " << name << " basis values..." << std::endl;
  
    const double width = 20.0;
    const double height = 10.0;
    
    typedef typename basis_t::cell_type cell_type;
    cfd::TriangularMeshBuilder meshBuilder(width, height, 2.0/15.0);
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
        assertZero(1.0 - sum);
      }
    }
  }

  template<typename basis_t>
  void testBasisDofs(const std::string& name)
  {
    std::cout << "Testing " << name << " basis function degrees-of-freedom..." << std::endl;
 
    const double width = 20.0;
    const double height = 10.0;
    
    typedef typename basis_t::cell_type cell_type;
    typedef typename cell_type::vertex_type vertex_type;
    typedef typename cfd::dof_map<cell_type> dof_map_t;
    typedef typename dof_map_t::local2global_map local2global_map;
    typedef typename local2global_map::key_type local_dof_t;
    typedef typename local2global_map::mapped_type global_dof_t;
    typedef typename std::map< global_dof_t, std::vector<local_dof_t> > global2local_map;

    cfd::TriangularMeshBuilder meshBuilder(width, height, 2.0/15.0);
    cfd::mesh<cell_type> m(meshBuilder.buildMesh());
    basis_t basis(m);

    cfd::dof_map_builder<cell_type> mapBuilder(m);
    mapBuilder.addFiniteElement(basis);

    dof_map_t dofMap(mapBuilder.getDofMap());
    std::map< global_dof_t, std::vector<local_dof_t> > global2local;

    for(typename dof_map_t::const_iterator dofIter(dofMap.begin()); dofIter!=dofMap.end(); ++dofIter)
    {
      global2local[dofIter->second].push_back(dofIter->first);
    }

    for(typename global2local_map::const_iterator dofIter(global2local.begin()); dofIter!=global2local.end(); ++dofIter)
    {
      if (dofIter->second.size() > 1)
      {
        const local_dof_t localDof(dofIter->second[0]);
        const vertex_type location = boost::get<0>(localDof)->getDofCoordinate(boost::get<1>(localDof), boost::get<2>(localDof));
        const double localDofValue = basis.evaluate_tensor(boost::get<1>(localDof), boost::get<2>(localDof), location).toScalar();

        for(unsigned dof = 0; dof < dofIter->second.size(); ++dof)
        {
          const local_dof_t coDof(dofIter->second[dof]);
          const vertex_type coDofLocation = boost::get<0>(coDof)->getDofCoordinate(boost::get<1>(coDof), boost::get<2>(coDof));
          const double coDofValue = basis.evaluate_tensor(boost::get<1>(coDof), boost::get<2>(coDof), coDofLocation).toScalar();
          assertTrue(location == coDofLocation);
          assertZero(localDofValue - coDofValue);
        }
      }
    }
  }

public:
  Tester();
  void run();
};
