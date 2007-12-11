#ifndef SIMPLE_CFD_DOF_MAP_HPP
#define SIMPLE_CFD_DOF_MAP_HPP

#include "simple_cfd_fwd.hpp"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <boost/lambda/lambda.hpp>

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

public:
  dof_map(const std::set<const finite_element_type*>& e, const local2global_map& m) : elements(e), mapping(m)
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
};

}

#endif
