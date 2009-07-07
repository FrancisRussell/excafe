#ifndef SIMPLE_CFD_BOUNDARY_CONDITION_HPP
#define SIMPLE_CFD_BOUBDARY_CONDITION_HPP

#include "simple_cfd_fwd.hpp"
#include "vertex.hpp"
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
  const Function<dimension, rank, double>& function;

public:
  BoundaryCondition(const SubDomain<dimension>& s, const Function<dimension, rank, double>& f) : subdomain(s), function(f)
  {
  }

  template<typename finite_element_t>
  void populateDirichletValues(FEVector<typename finite_element_t::cell_type>& boundaryValues, const finite_element_t& element)
  {
    BOOST_STATIC_ASSERT(rank == finite_element_t::rank);
    BOOST_STATIC_ASSERT(dimension == finite_element_t::dimension);

    typedef typename finite_element_t::cell_type cell_type;
    typedef dof_map<typename finite_element_t::cell_type> dof_map_type;
    typedef typename dof_map_type::dof_t dof_t;

    const dof_map_type dofMap(boundaryValues.getRowMappings());

    for(typename dof_map_type::const_iterator dofMapIter(dofMap.begin()); dofMapIter!=dofMap.end(); ++dofMapIter)
    {
      const dof_t dof(dofMapIter->first);
      const vertex<dimension> dofLocation(element.getDofCoordinate(boost::get<1>(dof), boost::get<2>(dof)));

      if (boost::get<0>(dof) == &element && subdomain.inside(dofLocation))
      {
        const Tensor<D, R, double> boundaryValue(function.evaluate(dofLocation));
        const Tensor<D, R, double> basisAtDofLocation(element.evaluate_tensor(boost::get<1>(dof),
                                                                              boost::get<2>(dof), dofLocation));
        const double dofValue = boundaryValue.colon_product(basisAtDofLocation).toScalar();
        boundaryValues.setValues(1, &dof, &dofValue);
      }
    }
    boundaryValues.assemble();
  }
};

}

#endif
