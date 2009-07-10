#ifndef SIMPLE_CFD_DOF_MAP_BUILDER_HPP
#define SIMPLE_CFD_DOF_MAP_BUILDER_HPP

#include <map>
#include <utility>
#include <set>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include "utility.hpp"
#include "simple_cfd_fwd.hpp"
#include "finite_element.hpp"

namespace cfd
{

template<typename C>
class dof_map_builder
{
private:
  typedef C cell_type;
  typedef finite_element<cell_type> finite_element_type;
  typedef std::map<boost::tuple<const finite_element_type*, cell_id, unsigned>, unsigned> local2global_map;

  const mesh<cell_type>& m;
  std::set<const finite_element_type*> elements;

public:
  dof_map_builder(const mesh<cell_type>& _m) : m(_m)
  {
  }

  void addFiniteElement(const finite_element<cell_type>& element)
  {
    elements.insert(&element);
  }

  dof_map<cell_type> getDofMap() const
  {
    const std::size_t dimension = m.getDimension();
    std::size_t counter = 0;
    local2global_map local2global;
    std::set<unsigned> boundaryDofs;

    // Iterate over all cells in order
    for(typename mesh<cell_type>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      // Iterate over finite elements
      for(typename std::set<const finite_element<cell_type>*>::const_iterator elementIter = elements.begin(); elementIter!=elements.end(); ++elementIter)
      {
        // Add in all degrees of freedom for current cell
        for(unsigned dof=0; dof<(*elementIter)->space_dimension(); ++dof)
        {
          local2global[boost::make_tuple(*elementIter, cellIter->getIndex(), dof)] = counter;
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
              const typename local2global_map::const_iterator globalDofIter = local2global.find(boost::make_tuple(*elementIter, *incidentIter, dofIter->second));
              assert(globalDofIter != local2global.end());
              local2global[boost::make_tuple(*elementIter, cellIter->getIndex(), dofIter->first)] = globalDofIter->second;
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

    // Record degrees of freedom that lie on a boundary
    const std::vector< std::pair<vertex_id, vertex_id> > boundary(m.getEdgeFacets());
    for(typename mesh<cell_type>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      for(typename std::set<const finite_element<cell_type>*>::const_iterator elementIter = elements.begin(); elementIter!=elements.end(); ++elementIter)
      {
        const std::vector<unsigned> localDofs((*elementIter)->getBoundaryDegreesOfFreedom(m, cellIter->getIndex(), boundary));
        for(std::vector<unsigned>::const_iterator localDofIter(localDofs.begin()); localDofIter!=localDofs.end(); ++localDofIter)
        {
          const typename local2global_map::const_iterator
            globalDofIter(local2global.find(boost::make_tuple(*elementIter, cellIter->getIndex(), *localDofIter)));
          assert(globalDofIter != local2global.end());
          boundaryDofs.insert(globalDofIter->second);
        }
      }
    }

    return dof_map<cell_type>(m, elements, local2global, boundaryDofs);
  }
};

}

#endif
