#ifndef SIMPLE_CFD_BOUNDARY_CONDITION_HPP
#define SIMPLE_CFD_BOUBDARY_CONDITION_HPP

#include "simple_cfd_fwd.hpp"
#include "vertex.hpp"
#include "dof.hpp"
#include "fe_vector.hpp"
#include "subdomain.hpp"
#include "function.hpp"
#include "cell_vertices.hpp"
#include <boost/static_assert.hpp>

namespace cfd
{

template<unsigned int D, unsigned int R>
class BoundaryCondition
{
public:
  static const unsigned int dimension = D;
  static const unsigned int rank = R;

private:
  const SubDomain<dimension>& subdomain;
  const Function<dimension>& function;

public:
  BoundaryCondition(const SubDomain<dimension>& s, const Function<dimension>& f) : subdomain(s), function(f)
  {
  }

  template<typename finite_element_t>
  void populateDirichletValues(FEVector<finite_element_t::dimension>& boundaryValues, const finite_element_t& element)
  {
    BOOST_STATIC_ASSERT(rank == finite_element_t::rank);
    BOOST_STATIC_ASSERT(dimension == finite_element_t::dimension);

    typedef typename finite_element_t::cell_type cell_type;
    typedef DofMap<finite_element_t::dimension> dof_map_type;
    typedef typename dof_map_type::dof_t dof_t;

    const dof_map_type dofMap(boundaryValues.getRowMappings());
    const Mesh<dimension>& m = dofMap.getMesh();

    for(typename dof_map_type::const_iterator dofMapIter(dofMap.begin()); dofMapIter!=dofMap.end(); ++dofMapIter)
    {
      const dof_t dof(dofMapIter->first);
      const vertex<dimension> dofLocationLocal(element.getDofCoordinateLocal(dof.getIndex()));
      const vertex<dimension> dofLocationGlobal(element.getDofCoordinateGlobal(m, dof.getCell(), dof.getIndex()));

      if (dof.getElement() == &element && subdomain.inside(dofLocationGlobal))
      {
        const CellVertices<dimension> vertices = m.getCoordinates(dof.getCell());
        const Tensor<D> boundaryValue(function.evaluate(dofLocationGlobal));
        const Tensor<D> basisAtDofLocation(element.evaluateTensor(vertices, dof.getIndex(), dofLocationLocal));
        const double dofValue = boundaryValue.colon_product(basisAtDofLocation);
        boundaryValues.setValues(1, &dof, &dofValue);
      }
    }
    boundaryValues.assemble();
  }
};

}

#endif
