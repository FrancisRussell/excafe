#ifndef SIMPLE_CFD_BOUNDARY_CONDITION_HPP
#define SIMPLE_CFD_BOUBDARY_CONDITION_HPP

#include "simple_cfd_fwd.hpp"
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
  FEVector<typename finite_element_t::cell_type> getDirichletValues(const dof_map<typename finite_element_t::cell_type>& dofMap, const finite_element_t& element)
  {
    BOOST_STATIC_ASSERT(rank == finite_element_t::rank);
    BOOST_STATIC_ASSERT(dimension == finite_element_t::dimension);

    typedef typename finite_element_t::cell_type cell_type;
    typedef dof_map<typename finite_element_t::cell_type> dof_map_type;

    dof_map<cell_type> subDomainDofMap(dofMap.extractSubDomainDofs(element, subdomain)); 

    FEVector<cell_type> dirichletValues(subDomainDofMap);

    for(typename dof_map_type::const_iterator dofMapIter(subDomainDofMap.begin()); dofMapIter!=subDomainDofMap.end(); ++dofMapIter)
    {
    }
  }
};

}

#endif
