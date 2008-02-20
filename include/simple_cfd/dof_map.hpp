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
private:
  typedef C cell_type;
  typedef finite_element<cell_type> finite_element_type;
  typedef std::map<boost::tuple<const finite_element_type*, cell_id, unsigned>, unsigned> local2global_map;
  const mesh<cell_type>& m;
  std::set<const finite_element_type*> elements;
  local2global_map mapping;
  std::set<unsigned> boundaryDofs;

public:
  dof_map()
  {
  }

  dof_map(const mesh<cell_type>& _m, const std::set<const finite_element_type*>& _elements, const local2global_map& _mapping, const std::set<unsigned> _boundaryDofs) : 
          m(_m), elements(_elements), mapping(_mapping), boundaryDofs(_boundaryDofs)
  {
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
    for(typename std::set<const finite_element_type*>::const_iterator elementIter(elements.begin()); elementIter!=elements.end(); ++elementIter)
      dofsPerCell += (*elementIter)->space_dimension();

    return dofsPerCell;
  }

  std::map< unsigned, std::set< boost::tuple<cell_id, unsigned> > > getBoundaryDegreesOfFreedom(const finite_element_type* const element) const
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

  unsigned getGlobalIndex(const boost::tuple<const finite_element_type*, cell_id, unsigned>& dof) const
  {
    const typename local2global_map::const_iterator mappingIter(mapping.find(dof));
    assert(mappingIter != mapping.end());
    return mappingIter->second;
  }

  SparsityPattern getSparsityPattern() const
  {
    const unsigned dofs = getDegreesOfFreedomCount();
    SparsityPattern pattern(dofs, dofs);
    const std::map<cell_id, cell_type> cells(m.getCells());

    for(typename std::map<cell_id, cell_type>::const_iterator cellIter(cells.begin()); cellIter!=cells.end(); ++cellIter)
      for(typename std::set<const finite_element_type*>::const_iterator testIter(elements.begin()); testIter!=elements.end(); ++testIter)
        for(typename std::set<const finite_element_type*>::const_iterator trialIter(elements.begin()); trialIter!=elements.end(); ++trialIter)
          for(unsigned test=0; test < (*testIter)->space_dimension(); ++test)
            for(unsigned trial=0; trial < (*trialIter)->space_dimension(); ++trial)
              pattern.insert(getGlobalIndex(boost::make_tuple(*testIter, cellIter->first, test)), getGlobalIndex(boost::make_tuple(*trialIter, cellIter->first, trial)));

    return pattern;
  }
};

}

#endif
