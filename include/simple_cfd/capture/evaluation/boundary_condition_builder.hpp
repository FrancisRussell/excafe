#ifndef SIMPLE_CFD_CAPTURE_BOUNDARY_CONDITION_BUILDER_HPP
#define SIMPLE_CFD_CAPTURE_BOUNDARY_CONDITION_BUILDER_HPP

#include <cstddef>
#include <set>
#include <simple_cfd/simple_cfd_fwd.hpp>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/dof.hpp>
#include <simple_cfd/dof_map.hpp>
#include <simple_cfd/boundary_condition_list.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class BoundaryConditionBuilder
{
private:
  static const std::size_t dimension = D;
  typedef Dof<dimension> dof_t;
  
  Mesh<dimension>& mesh;

  std::set<dof_t> getDofs(const FiniteElement<dimension>& element, const BoundaryConditionList<dimension>& condition) const
  {
    std::set<dof_t> dofs;

    // We iterate over each facet on each cell  instead of iterating over all facets directly. This allows
    // us to handle degrees-of-freedom shared between cells.
    for(typename Mesh<dimension>::global_iterator cIter(mesh.global_begin(dimension)); cIter!=mesh.global_end(dimension); ++cIter)
    {
      for(typename Mesh<dimension>::local_iterator fIter(mesh.local_begin(*cIter, dimension-1)); 
        fIter!=mesh.local_end(*cIter, dimension-1); ++fIter)
      {
        if (condition.applies(mesh.getFacetLabel(*fIter)))
        {
          const std::set<dof_t> dofsOnFacet = element.getDofsOnEntity(mesh.getTopology(), cIter->getIndex(), *fIter);
          dofs.insert(dofsOnFacet.begin(), dofsOnFacet.end());
        }
      }
    }
    return dofs;
  }

  DiscreteField<dimension> getBoundaryField(const DofMap<dimension>& map, const BoundaryConditionList<dimension>& condition) const
  {
    const FiniteElement<dimension>& element = *map.getFiniteElement();
    DiscreteField<dimension> boundaryField(map);
    std::map<std::size_t, int> priorityMap;


    for(typename Mesh<dimension>::global_iterator cIter(mesh.global_begin(dimension)); cIter!=mesh.global_end(dimension); ++cIter)
    {
      for(typename Mesh<dimension>::local_iterator fIter(mesh.local_begin(*cIter, dimension-1)); 
        fIter!=mesh.local_end(*cIter, dimension-1); ++fIter)
      {
        const int label= mesh.getFacetLabel(*fIter);
        if (condition.applies(label))
        {
          const std::set<dof_t> dofsOnFacet = element.getDofsOnEntity(mesh.getTopology(), cIter->getIndex(), *fIter);

          for(typename std::set<dof_t>::const_iterator dofIter(dofsOnFacet.begin()); dofIter!=dofsOnFacet.end(); ++dofIter)
          {
            const std::size_t index = map.getGlobalIndex(*dofIter);
            const std::map<std::size_t, int>::const_iterator priorityIter = priorityMap.find(index);

            if (priorityIter == priorityMap.end() || priorityIter->second < condition.getPriority(label))
            {
              priorityMap[index] = condition.getPriority(label);

              const vertex<dimension> dofLocation = element.getDofCoordinateGlobal(mesh, dofIter->getCell(), dofIter->getIndex());
              const Tensor<dimension> boundaryValue = condition.getValue(dofLocation, label);
              const std::size_t tensorIndex = element.getTensorIndex(mesh, dofIter->getCell(), dofIter->getIndex());

              //FIXME: the tensor index operator only works on rank 1 fields
              const double dofCoeff = boundaryValue(tensorIndex);

              boundaryField.setValues(1, &(*dofIter), &dofCoeff);
            }
          }
        }
      }
    }

    return boundaryField;
  }

public:
  BoundaryConditionBuilder(Mesh<dimension>& _mesh) : mesh(_mesh)
  {
  }

  DiscreteField<dimension> getBoundaryValues(const DofMap<dimension>& space, const BoundaryConditionList<dimension>& condition)
  {
    const std::set<dof_t> dofs = getDofs(*space.getFiniteElement(), condition);
    const DofMap<dimension> boundaryDofMap(space.intersect(dofs));
    return getBoundaryField(boundaryDofMap, condition);
  }
};

}

}

#endif
