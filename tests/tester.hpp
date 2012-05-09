#include <excafe/mesh.hpp>
#include <excafe/triangular_mesh_builder.hpp>
#include <excafe/dof_map_builder.hpp>
#include <excafe/mesh_entity.hpp>
#include <excafe/quadrature_points.hpp>
#include <excafe/cell_vertices.hpp>
#include <excafe/dof.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <boost/array.hpp>

class Tester
{
private:
  const double epsilon;
  typedef excafe::TriangularCell cell_type;
  typedef cell_type::vertex_type vertex_type;

  void assertTrue(const bool b);
  void assertFalse(const bool b);
  void assertEqual(const double d, const double d2);

  void testQuadrature();  
  void testQuadrature(const std::map<double, double>& q);  
  void testTriangleQuadrature();  
  void testLinearBasis();
  void testQuadraticBasis();
  void testLinearBasisDofs();
  void testQuadraticBasisDofs();

  template<typename basis_t>
  void testBasis(const basis_t& basis, const std::string& name)
  {
    std::cout << "Testing " << name << " basis values..." << std::endl;
  
    const double width = 20.0;
    const double height = 10.0;
    
    typedef typename basis_t::cell_type cell_type;
    excafe::TriangularMeshBuilder meshBuilder(width, height, 2.0/15.0);
    excafe::Mesh<cell_type::dimension> m(meshBuilder.buildMesh());
    const std::size_t dimension = m.getDimension();

    const std::size_t degree = 5;
    boost::array<std::size_t, 2> degrees;
    std::fill(degrees.begin(), degrees.end(), degree);

    excafe::QuadraturePoints<2> quadrature(m.getReferenceCell()->getQuadrature(degrees));
    const excafe::MeshEntity localCell(dimension, 0);
  
    for(typename excafe::Mesh<cell_type::dimension>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      const int dofs = basis.spaceDimension();

      for(excafe::QuadraturePoints<2>::const_iterator wIter(quadrature.begin(localCell)); wIter!=quadrature.end(localCell); ++wIter)
      {
        double sum = 0.0;
        for(int i=0; i<dofs; ++i)
          sum += basis.evaluateTensor(i, wIter->first);

        assertEqual(1.0, sum);
      }
    }
  }

  template<typename basis_t>
  void testBasisDofs(const basis_t& basis, const std::string& name)
  {
    std::cout << "Testing " << name << " basis function degrees-of-freedom..." << std::endl;
 
    const double width = 20.0;
    const double height = 10.0;
    
    typedef typename basis_t::cell_type cell_type;
    typedef typename cell_type::vertex_type vertex_type;
    typedef typename excafe::DofMap<cell_type::dimension> dof_map_t;
    typedef typename dof_map_t::local2global_map local2global_map;
    typedef typename local2global_map::key_type local_dof_t;
    typedef typename local2global_map::mapped_type global_dof_t;
    typedef typename std::map< global_dof_t, std::vector<local_dof_t> > global2local_map;

    excafe::TriangularMeshBuilder meshBuilder(width, height, 2.0/15.0);
    excafe::Mesh<cell_type::dimension> m(meshBuilder.buildMesh());

    excafe::DofMapBuilder<cell_type::dimension> mapBuilder(m);
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
        const vertex_type location = localDof.getElement()->getDofCoordinateLocal(localDof.getIndex());
        const double localDofValue = basis.evaluateTensor(localDof.getIndex(), location);

        for(unsigned dof = 0; dof < dofIter->second.size(); ++dof)
        {
          const local_dof_t coDof(dofIter->second[dof]);
          const vertex_type coDofLocation = coDof.getElement()->getDofCoordinateLocal(coDof.getIndex());
          const double coDofValue = basis.evaluateTensor(coDof.getIndex(), coDofLocation);
          assertTrue(m.referenceToPhysical(localDof.getCell(), location) == m.referenceToPhysical(coDof.getCell(), coDofLocation));
          assertEqual(localDofValue, coDofValue);
        }
      }
    }
  }

public:
  Tester();
  void run();
};
