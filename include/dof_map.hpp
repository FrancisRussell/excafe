#ifndef SIMPLE_CFD_DOF_MAP_HPP
#define SIMPLE_CFD_DOF_MAP_HPP

#include "simple_cfd_fwd.hpp"
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
private:
  typedef C cell_type;
  typedef finite_element<cell_type> finite_element_type;
  typedef std::map<boost::tuple<const finite_element_type*, cell_id, unsigned>, unsigned> local2global_map;
  std::set<const finite_element_type*> elements;
  local2global_map mapping;
  std::set<unsigned> boundaryDofs;

public:
  dof_map()
  {
  }

  dof_map(const std::set<const finite_element_type*>& _elements, const local2global_map& _mapping, const std::set<unsigned> _boundaryDofs) : 
          elements(_elements), mapping(_mapping), boundaryDofs(_boundaryDofs)
  {
  }

  std::size_t getMappingSize() const
  {
    return mapping.size();
  }

  std::size_t getDegreesOfFreedom() const
  {
    std::set<unsigned> dofs;
    for(typename local2global_map::const_iterator mapIter(mapping.begin()); mapIter!=mapping.end(); ++mapIter)
      dofs.insert(mapIter->second);
    return dofs.size();
  }

  std::size_t getBoundaryDegreesOfFreedom() const
  {
    return boundaryDofs.size();
  }

  unsigned getGlobalIndex(const boost::tuple<const finite_element_type*, cell_id, unsigned>& dof) const
  {
    const typename local2global_map::const_iterator mappingIter(mapping.find(dof));
    assert(mappingIter != mapping.end());
    return mappingIter->second;
  }
};

}

#endif
