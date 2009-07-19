#ifndef SIMPLE_CFD_DOF_MAP_BUILDER_HPP
#define SIMPLE_CFD_DOF_MAP_BUILDER_HPP

#include <map>
#include <utility>
#include <set>
#include <algorithm>
#include <cassert>
#include <boost/foreach.hpp>
#include "simple_cfd_fwd.hpp"
#include "finite_element.hpp"
#include "dof_map.hpp"
#include "dof.hpp"

namespace cfd
{

template<typename C>
class DofMapBuilder
{
private:
  typedef C cell_type;
  static const std::size_t dimension = cell_type::dimension;
  typedef FiniteElement<dimension> finite_element_t;
  typedef Dof<dimension> dof_t;
  typedef std::map<dof_t, unsigned> local2global_map;

  const Mesh<dimension>& m;
  std::set<const finite_element_t*> elements;

public:
  DofMapBuilder(const Mesh<dimension>& _m) : m(_m)
  {
  }

  void addFiniteElement(const finite_element_t& element)
  {
    elements.insert(&element);
  }

  DofMap<cell_type> getDofMap() const
  {
    std::size_t counter = 0;
    local2global_map local2global;

    for(std::size_t d=0; d<=dimension; ++d)
    {
      for(typename Mesh<dimension>::global_iterator eIter(m.global_begin(d)); eIter!=m.global_end(d); ++eIter)
      {
        for(typename std::set<const finite_element_t*>::const_iterator elementIter = elements.begin(); elementIter!=elements.end(); ++elementIter)
        {
          std::set<dof_t> dofsOnEntity;

          if (d < dimension)
          {
            for(typename Mesh<dimension>::local_iterator cellIter(m.local_begin(*eIter, dimension)); cellIter!=m.local_end(*eIter, dimension); ++cellIter)
            {
              const std::set<dof_t> cellDofsOnEntity = (*elementIter)->getDofsOnEntity(m.getTopology(), cellIter->getIndex(), *eIter);
              dofsOnEntity.insert(cellDofsOnEntity.begin(), cellDofsOnEntity.end());
            }
          }
          else
          {
            // When we're iterating over cells, we must avoid looking at other cells
            dofsOnEntity = (*elementIter)->getDofsOnEntity(m.getTopology(), eIter->getIndex(), *eIter);
          }

          const std::vector< std::set<dof_t> > identicalDofsList = (*elementIter)->resolveIdenticalDofs(m, *eIter, dofsOnEntity);
          BOOST_FOREACH(const std::set<dof_t>& identicalDofs, identicalDofsList)
          {
            assert(!identicalDofs.empty());
            BOOST_FOREACH(const dof_t& dof, identicalDofs)
            {
              local2global[dof] = counter;
            }
            ++counter;
          }
        }
      }
    }

    return DofMap<cell_type>(m, elements, local2global);
  } 
};

}

#endif
