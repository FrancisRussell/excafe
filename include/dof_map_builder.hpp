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
  unsigned counter;
  std::set<const finite_element_type*> elements;
  local2global_map local2global;

  void addDof(const boost::tuple<finite_element<cell_type>*, cell_id, unsigned>& dof)
  {
  }

public:
  dof_map_builder(const mesh<cell_type>& _m) : m(_m), counter(0)
  {
  }

  void addFiniteElement(const finite_element<cell_type>& element)
  {
    elements.insert(&element);
  }

  dof_map<cell_type> getDofMap() const
  {
    return dof_map<cell_type>(elements, local2global);
  }

  void handleCells(const std::map<cell_id, cell_type>& cells)
  {
/*  // Find cells which have vertices in common
    std::map< vertex_id, std::set<cell_id> > adjacencyByVertex;

    for(typename std::map<cell_id, cell_type>::const_iterator cellIter = cells.begin(); cellIter!=cells.end(); ++cellIter)
    {
      const std::vector<vertex_id> indices(cellIter->second.getIndices());
      for(std::vector<vertex_id>::const_iterator vertexIter = indices.begin(); vertexIter != indices.end(); ++vertexIter)
      {
        adjacencyByVertex[*vertexIter].insert(cellIter->first);
      }
    }

    // Create a cell adjacency set
    typedef std::set< std::pair<cell_id, cell_id>, unordered_pair_compare<cell_id> > adjacency_set;
    adjacency_set cellAdjacency;
    for(std::map< vertex_id, std::set<cell_id> >::const_iterator abvIter = adjacencyByVertex.begin(); abvIter!=adjacencyByVertex.end(); ++abvIter)
    {
      const std::set<cell_id> adjacentCells(abvIter->second);
      for(std::set<cell_id>::iterator cellIter = adjacentCells.begin(); cellIter!=adjacentCells.end(); ++cellIter)
      {
        for(std::set<cell_id>::iterator cellIter2 = adjacentCells.begin(); cellIter2!=adjacentCells.end(); ++cellIter2)
          cellAdjacency.insert(std::make_pair(*cellIter, *cellIter2));
      }
    }
*/
    // Iterate over all cells in order
    for(typename std::map<cell_id, cell_type>::const_iterator cellIter = cells.begin(); cellIter!=cells.end(); ++cellIter)
    {
      // Iterate over finite elements
      for(typename std::set<const finite_element<cell_type>*>::const_iterator elementIter = elements.begin(); elementIter!=elements.end(); ++elementIter)
      {
        // Add in all degrees of freedom for current cell
        for(unsigned dof=0; dof<(*elementIter)->space_dimension(); ++dof)
        {
          local2global[boost::make_tuple(*elementIter, cellIter->first, dof)] = counter;
          ++counter;
        }

        // Iterate over cells incident to current one
        const std::set<cell_id> incident(m.getCellIncidentCells(cellIter->first));
        for(std::set<cell_id>::const_iterator incidentIter(incident.begin()); incidentIter != incident.end(); ++incidentIter)
        {
          const std::vector< std::pair<unsigned, unsigned> > commonDofs = (*elementIter)->getCommonDegreesOfFreedom(cellIter->first, *incidentIter);

          // Iterate over map of degrees of freedom of current cell to incident cell
          for(std::vector< std::pair<unsigned, unsigned> >::const_iterator dofIter(commonDofs.begin()); dofIter!=commonDofs.end(); ++dofIter)
          {
            // If we are incident with a lower numbered cell, we replace our mapping with the one for that
            // cell 
            if (*incidentIter < cellIter->first)
            {
              const typename local2global_map::const_iterator globalDofIter = local2global.find(boost::make_tuple(*elementIter, *incidentIter, dofIter->second));
              assert(globalDofIter != local2global.end());
              local2global[boost::make_tuple(*elementIter, cellIter->first, dofIter->first)] = globalDofIter->second;
            }
          }
        }
      }
    }

    // Now we must renumber the global degrees of freedom so they aren't sparse
    unsigned newCounter = 0;
    std::map<unsigned, unsigned> newMapping;
    for(typename local2global_map::iterator mapIter = local2global.begin(); mapIter!=local2global.end(); ++mapIter)
    {
      const std::map<unsigned, unsigned>::const_iterator newMappingIter = newMapping.find(mapIter->second);
      if (newMappingIter == newMapping.end())
      {
        newMapping[mapIter->second] = newCounter;
        mapIter->second = newCounter;
        ++newCounter;
      }
      else
      {
        mapIter->second = newMappingIter->second;
      }
    }
  }
};

}

#endif
