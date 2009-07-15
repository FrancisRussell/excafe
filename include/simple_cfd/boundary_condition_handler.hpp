#ifndef SIMPLE_CFD_BOUNDARY_CONDITION_HANDLER_HPP
#define SIMPLE_CFD_BOUNDARY_CONDITION_HANDLER_HPP

#include <cstddef>
#include "simple_cfd_fwd.hpp"
#include "boundary_condition2.hpp"
#include "dof.hpp"

namespace cfd
{

class BoundaryConditionHandler
{
private:
  static const std::size_t dimension = 2;
  typedef Mesh<dimension> mesh_t;
  typedef Dof<dimension> local_dof_t;
  mesh_t& m;

  std::size_t findCell(const MeshEntity& entity) const
  {
    if (entity.getDimension() == m.getDimension())
      return entity.getIndex();

    std::vector<std::size_t> cellIndices(m.getTopology().getIndices(entity, m.getDimension()));
    assert(cellIndices.size() > 0);
    return cellIndices.front();
  }

  template<typename finite_element_t>
  std::set<local_dof_t> getEntityDofs(const finite_element_t& element, const MeshEntity& entity) const
  {
    const std::size_t cellId = findCell(entity);

    // Find all degrees-of-freedom on the entity
    std::set<local_dof_t> localDofs(element.getDegreesOfFreedom(m.getTopology(), cellId, entity));

    // We also need to look at degrees-of-freedom on the entity tied to entities of lower topological
    // dimension
    for(std::size_t dofEntityDim=0; dofEntityDim<entity.getDimension(); ++dofEntityDim)
    {
      std::vector<std::size_t> entityIndices(m.getIndices(entity, dofEntityDim));

      // Iterate over all mesh entities of dimension dofEntityDim of the facet
      for(std::vector<std::size_t>::const_iterator dofEntityIter(entityIndices.begin()); dofEntityIter!=entityIndices.end(); ++dofEntityIter)
      {
        const std::set<local_dof_t> lowerDofs(element.getDegreesOfFreedom(m.getTopology(), cellId, MeshEntity(dofEntityDim, *dofEntityIter)));
        localDofs.insert(lowerDofs.begin(), lowerDofs.end());
      }
    }

    return localDofs;
  }

public:
  BoundaryConditionHandler(mesh_t& _m) : m(_m)
  {
  }

  template<typename finite_element_t>
  void handleBoundaryCondition(const finite_element_t& element, 
    const BoundaryCondition2<finite_element_t::dimension, finite_element_t::rank>& condition)
  {
    const std::size_t dimension = m.getDimension();
    std::map< std::size_t, std::set<local_dof_t> > facetDofMap; 

    for(mesh_t::global_iterator facetIter(m.global_begin(dimension-1)); facetIter!=m.global_end(dimension-1); ++facetIter)
    {
      const int label = m.getFacetLabel(*facetIter);
      if (condition.applies(m.getGeometry(), *facetIter, label))
      {
        // Find all degrees-of-freedom on the facet and of lower-dimension entities on the facet
        std::set<local_dof_t> facetDofs(getEntityDofs(element, *facetIter));

        // Filter degrees-of-freedom if they only apply to one tensor index
        if (condition.constrainsIndex(m.getTopology(), *facetIter, label))
        {
          const std::size_t constrainedIndex = condition.constrainedIndex(m.getTopology(), *facetIter, label);
          std::set<local_dof_t> filteredDofs;

          for(std::set<local_dof_t>::const_iterator dofIter(facetDofs.begin()); dofIter!=facetDofs.end(); ++dofIter)
          {
            if (element.getTensorIndex(dofIter->getCell(), dofIter->getIndex()) == constrainedIndex)
              filteredDofs.insert(*dofIter);
          }

          facetDofMap[facetIter->getIndex()] = filteredDofs;
        }
        else
        {
          facetDofMap[facetIter->getIndex()] = facetDofs;
        }
      }

      // Now we've collected all the dofs, we need to integrate them over each facet
    }
  }
};

}

#endif
