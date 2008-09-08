#ifndef SIMPLE_CFD_DOF_MAP_HPP
#define SIMPLE_CFD_DOF_MAP_HPP

#include "simple_cfd_fwd.hpp"
#include "numeric/sparsity_pattern.hpp"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <cassert>
#include <map>
#include <set>

namespace cfd
{

template<typename C>
class dof_map
{
public:
  typedef C cell_type;
  typedef finite_element<cell_type> finite_element_t;
  typedef boost::tuple<const finite_element_t*, cell_id, unsigned>  dof_t;
  typedef std::map<dof_t, unsigned> local2global_map;

private:
  const mesh<cell_type>& m;
  std::set<const finite_element_t*> elements;
  local2global_map mapping;
  std::set<unsigned> boundaryDofs;

public:
  dof_map()
  {
  }

  dof_map(const mesh<cell_type>& _m, const std::set<const finite_element_t*>& _elements, const local2global_map& _mapping, const std::set<unsigned> _boundaryDofs) : 
          m(_m), elements(_elements), mapping(_mapping), boundaryDofs(_boundaryDofs)
  {
  }

  const mesh<cell_type>& getMesh() const
  {
    return m;
  }

  std::set<const finite_element_t*> getFiniteElements() const
  {
    return elements;
  }

  std::size_t getMappingSize() const
  {
    return mapping.size();
  }

  std::size_t getDegreesOfFreedomCount() const
  {
    std::set<unsigned> dofs;
    for(typename local2global_map::const_iterator mapIter(mapping.begin()); mapIter!=mapping.end(); ++mapIter)
      dofs.insert(mapIter->second);
    return dofs.size();
  }

  std::size_t getBoundaryDegreesOfFreedomCount() const
  {
    return boundaryDofs.size();
  }

  std::size_t getDofsPerCell() const
  {
    std::size_t dofsPerCell = 0;
    for(typename std::set<const finite_element_t*>::const_iterator elementIter(elements.begin()); elementIter!=elements.end(); ++elementIter)
      dofsPerCell += (*elementIter)->space_dimension();

    return dofsPerCell;
  }

  std::map< unsigned, std::set< boost::tuple<cell_id, unsigned> > > getBoundaryDegreesOfFreedom(const finite_element_t* const element) const
  {
    std::map< unsigned, std::set<boost::tuple<cell_id, unsigned> > > global2local;
    for(typename local2global_map::const_iterator mappingIter(mapping.begin()); mappingIter != mapping.end(); ++mappingIter)
    {
      if (boundaryDofs.find(mappingIter->second) != boundaryDofs.end())
      {
        // TODO: Work out why tuple.get<N>() won't compile
        if (boost::get<0>(mappingIter->first) == element)
          global2local[mappingIter->second].insert(boost::make_tuple(boost::get<1>(mappingIter->first), boost::get<2>(mappingIter->first)));
      }
    }
    return global2local;
  }

  dof_map extractDofs(const finite_element_t* element) const
  {
    std::set<const finite_element_t*> newElements;
    std::set<unsigned> newDofs;
    local2global_map newMapping;

    // We only define a mapping for a single element
    newElements.insert(element);

    // Copy across relevant mappings and get all global dofs corresponding to element so we can filter
    // boundary dofs;
    for(typename local2global_map::const_iterator mappingIter=mapping.begin(); mappingIter!=mapping.end(); ++mappingIter)
    {
      if (boost::get<0>(mappingIter->first) == element)
      {
        newDofs.insert(mappingIter->second);
        newMapping.insert(*mappingIter);
      }
    }

    // Work out the subset of dofs on the boundary
    std::set<unsigned> newBoundaryDofs;
    std::set_intersection(boundaryDofs.begin(), boundaryDofs.end(), newDofs.begin(), newDofs.end(), newBoundaryDofs.begin());

    return dof_map(m, newElements, newMapping, newBoundaryDofs);
  }

  dof_map contiguous() const
  {
    // Create set of current dof values
    std::set<unsigned> currentDofs;
    for(typename local2global_map::const_iterator mappingIter=mapping.begin(); mappingIter!=mapping.end(); ++mappingIter)
      currentDofs.insert(mappingIter->second);

    // Create the mapping between old and new dofs
    unsigned newDof = 0;
    std::map<unsigned, unsigned> remapping;
    for(std::set<unsigned>::const_iterator currentDofIter(currentDofs.begin()); currentDofIter!=currentDofs.end(); ++currentDofIter)
      remapping[*currentDofIter] = newDof++;

    // Remap the mappings
    local2global_map newMapping;
    for(typename local2global_map::const_iterator mappingIter=mapping.begin(); mappingIter!=mapping.end(); ++mappingIter)
      newMapping[mappingIter->first] = remapping[mappingIter->second];

    // Remap the boundary dofs
    std::set<unsigned> newBoundaryDofs;
    for(std::set<unsigned>::const_iterator boundaryDofIter(boundaryDofs.begin()); boundaryDofIter != boundaryDofs.end(); ++boundaryDofIter)
      newBoundaryDofs.insert(remapping[*boundaryDofIter]);

    return dof_map(m, elements, newMapping, newBoundaryDofs);
  }

  unsigned getGlobalIndex(const dof_t& dof) const
  {
    const typename local2global_map::const_iterator mappingIter(mapping.find(dof));
    assert(mappingIter != mapping.end());
    return mappingIter->second;
  }
};

}

#endif
