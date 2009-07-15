#ifndef SIMPLE_CFD_DOF_MAP_BUILDER_HPP
#define SIMPLE_CFD_DOF_MAP_BUILDER_HPP

#include <map>
#include <utility>
#include <set>
#include <algorithm>
#include "utility.hpp"
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
    const std::size_t mesh_dimension = m.getDimension();
    std::size_t counter = 0;
    local2global_map local2global;

    // Iterate over all cells in order
    for(typename Mesh<dimension>::global_iterator cellIter(m.global_begin(mesh_dimension)); cellIter!=m.global_end(mesh_dimension); ++cellIter)
    {
      // Iterate over finite elements
      for(typename std::set<const finite_element_t*>::const_iterator elementIter = elements.begin(); elementIter!=elements.end(); ++elementIter)
      {
        // Add in all degrees of freedom for current cell
        for(unsigned dof=0; dof<(*elementIter)->spaceDimension(); ++dof)
        {
          local2global[dof_t(*elementIter, cellIter->getIndex(), dof)] = counter;
          ++counter;
        }

        // Iterate over cells incident to current one
        const std::set<cell_id> incident(m.getCellIncidentCells(cellIter->getIndex()));
        for(std::set<cell_id>::const_iterator incidentIter(incident.begin()); incidentIter != incident.end(); ++incidentIter)
        {
          const std::vector< std::pair<unsigned, unsigned> > commonDofs =
            (*elementIter)->getCommonDegreesOfFreedom(m, cellIter->getIndex(), *incidentIter);

          // Iterate over map of degrees of freedom of current cell to incident cell
          for(std::vector< std::pair<unsigned, unsigned> >::const_iterator dofIter(commonDofs.begin()); dofIter!=commonDofs.end(); ++dofIter)
          {
            // If we are incident with a lower numbered cell, we replace our mapping with the one for that
            // cell 
            if (*incidentIter < cellIter->getIndex())
            {
              const typename local2global_map::const_iterator globalDofIter = local2global.find(dof_t(*elementIter, *incidentIter, dofIter->second));
              assert(globalDofIter != local2global.end());
              local2global[dof_t(*elementIter, cellIter->getIndex(), dofIter->first)] = globalDofIter->second;
            }
          }
        }
      }
    }

    // Now we must renumber the global degrees of freedom so they aren't sparse
    // We do this in three phases because we want to try to preserve the order of the dofs,
    // keeping the matrix more diagonally dominant (maybe).
    // TODO: Implement some sort of renumbering scheme

    std::set<unsigned> dofs;
    for(typename local2global_map::const_iterator mapIter = local2global.begin(); mapIter!=local2global.end(); ++mapIter)
      dofs.insert(mapIter->second);

    unsigned newCounter=0;
    std::map<unsigned, unsigned> renumbering;
    for(std::set<unsigned>::const_iterator dofIter(dofs.begin()); dofIter!=dofs.end(); ++dofIter)
    {
      renumbering[*dofIter] = newCounter;
      ++newCounter;
    }

    for(typename local2global_map::iterator mapIter = local2global.begin(); mapIter!=local2global.end(); ++mapIter)
    {
      const std::map<unsigned, unsigned>::const_iterator renumberingIter = renumbering.find(mapIter->second);
      assert(renumberingIter != renumbering.end());
      mapIter->second = renumberingIter->second;
    }

    return DofMap<cell_type>(m, elements, local2global);
  }
};

}

#endif
