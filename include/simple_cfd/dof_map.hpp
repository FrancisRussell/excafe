#ifndef SIMPLE_CFD_DOF_MAP_HPP
#define SIMPLE_CFD_DOF_MAP_HPP

#include <cassert>
#include <map>
#include <set>
#include <iterator>
#include <utility>
#include <cstddef>
#include "simple_cfd_fwd.hpp"
#include "exception.hpp"
#include "dof.hpp"
#include "numeric/sparsity_pattern.hpp"

namespace cfd
{

template<std::size_t D>
class DofMap
{
public:
  static const std::size_t dimension = D;
  typedef FiniteElement<dimension> finite_element_t;
  typedef Dof<dimension> dof_t;
  typedef std::map<dof_t, unsigned> local2global_map;

  typedef typename local2global_map::const_iterator const_iterator;

private:
  const Mesh<dimension>* m;
  std::set<const finite_element_t*> elements;
  local2global_map mapping;

  void makeContiguous()
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
    for(typename local2global_map::iterator mappingIter=mapping.begin(); mappingIter!=mapping.end(); ++mappingIter)
      mappingIter->second = remapping[mappingIter->second];
  }

public:
  DofMap()
  {
  }

  DofMap(const DofMap& d) : m(d.m), elements(d.elements), mapping(d.mapping)
  {
  }

  DofMap(const Mesh<dimension>& _m, const std::set<const finite_element_t*>& _elements, const local2global_map& _mapping) : 
          m(&_m), elements(_elements), mapping(_mapping)
  {
  }
  
  bool isComposite() const
  {
    return elements.size() > 1;
  }

  const finite_element_t* getFiniteElement() const
  {
    assert(!isComposite());
    return *elements.begin();
  }

  const_iterator begin() const
  {
    return mapping.begin();
  }

  const_iterator end() const
  {
    return mapping.end();
  }

  bool operator+=(const DofMap& map)
  {
    if (m != map.m)
    {
      CFD_EXCEPTION("Attempted to add two DofMaps defined on different meshes");
    }
    //TODO: implement me!
  }

  bool operator==(const DofMap& map) const
  {
    return m == map.m &&
    elements == map.elements &&
    mapping == map.mapping;
  }

  const Mesh<dimension>& getMesh() const
  {
    return *m;
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

  std::size_t getDofsPerCell() const
  {
    std::size_t dofsPerCell = 0;
    for(typename std::set<const finite_element_t*>::const_iterator elementIter(elements.begin()); elementIter!=elements.end(); ++elementIter)
      dofsPerCell += (*elementIter)->space_dimension();

    return dofsPerCell;
  }

  DofMap extractDofs(const finite_element_t* element) const
  {
    std::set<const finite_element_t*> newElements;
    local2global_map newMapping;

    // We only define a mapping for a single element
    newElements.insert(element);

    // Copy across relevant mappings and get all global dofs corresponding to element so we can filter
    // boundary dofs;
    for(typename local2global_map::const_iterator mappingIter=mapping.begin(); mappingIter!=mapping.end(); ++mappingIter)
    {
      if (mappingIter->first.getElement() == element)
      {
        newMapping.insert(*mappingIter);
      }
    }

    DofMap result(*m, newElements, newMapping);
    result.makeContiguous();
    return result;
  }

  std::pair<DofMap, DofMap> splitHomogeneousDirichlet(const std::vector< std::pair<const finite_element_t*, const SubDomain<dimension>*> >& dirichletConditions)
  {
    local2global_map homogeneous;
    std::set<const finite_element_t*> homogeneousElements;

    local2global_map dirichlet;
    std::set<const finite_element_t*> dirichletElements;

    for(typename local2global_map::const_iterator mappingIter=mapping.begin(); mappingIter!=mapping.end(); ++mappingIter)
    {
      const dof_t dof(mappingIter->first);
      bool isDirichlet = false;

      for(typename std::vector< std::pair<const finite_element_t*, const SubDomain<dimension>*> >::const_iterator dirichletIter = dirichletConditions.begin(); 
          dirichletIter != dirichletConditions.end(); ++dirichletIter)
      {
        if (dirichletIter->first == dof.getElement() &&
          dirichletIter->second->inside(dof.getElement()->getDofCoordinateGlobal(*m, dof.getCell(), dof.getIndex())))
        {
          dirichlet.insert(*mappingIter);
          dirichletElements.insert(dof.getElement());

          isDirichlet = true;
        }
      }

      if (!isDirichlet) 
      {
        homogeneous.insert(*mappingIter);
        homogeneousElements.insert(dof.getElement());
      }
    }

    DofMap homogeneousMap(*m, homogeneousElements, homogeneous);
    DofMap dirichletMap(*m, dirichletElements, dirichlet);

    homogeneousMap.makeContiguous();
    dirichletMap.makeContiguous();

    return std::make_pair(homogeneousMap, dirichletMap);
  }

  std::vector<int> getIndices(const DofMap& map) const
  {
    std::vector<int> result(getDegreesOfFreedomCount());
    
    for(typename local2global_map::const_iterator mappingIter=mapping.begin(); mappingIter!=mapping.end(); ++mappingIter)
      result[mappingIter->second] = map.getGlobalIndex(mappingIter->first);

    return result;
  }

  unsigned getGlobalIndex(const dof_t& dof) const
  {
    const typename local2global_map::const_iterator mappingIter(mapping.find(dof));
    assert(mappingIter != mapping.end());
    return mappingIter->second;
  }

  int getGlobalIndexWithMissingAsNegative(const dof_t& dof) const
  {
    const typename local2global_map::const_iterator mappingIter(mapping.find(dof));
    if (mappingIter != mapping.end())
    {
      return mappingIter->second;
    }
    else
    {
      return -1;
    }
  }
};

}

#endif
