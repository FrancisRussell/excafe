#include <simple_cfd/cse/sop_map.hpp>
#include <simple_cfd/cse/polynomial_index.hpp>
#include <vector>
#include <utility>
#include <cstddef>
#include <map>
#include <cassert>

namespace cfd
{

namespace cse
{

SOPMap::const_iterator SOPMap::begin() const
{
  return map.begin();
}

SOPMap::const_iterator SOPMap::end() const
{
  return map.end();
}

std::vector<PolynomialIndex> SOPMap::reserveIndices(const std::size_t count)
{
  std::vector<PolynomialIndex> result;
  result.reserve(count);

  for(std::size_t i=0; i<count; ++i)
  {
    const PolynomialIndex index = newIndex();
    result.push_back(index);
    map.insert(std::make_pair(index, SOP()));
  }
  return result;
}

PolynomialIndex SOPMap::addSOP(const SOP& sop)
{
  const PolynomialIndex index = newIndex();
  map.insert(std::make_pair(index, sop));
  return index;
}

SOP& SOPMap::operator[](const PolynomialIndex& i)
{
  const std::map<PolynomialIndex, SOP>::iterator iter = map.find(i);
  assert(iter != map.end());
  return iter->second;
}

const SOP& SOPMap::operator[](const PolynomialIndex& i) const
{
  const std::map<PolynomialIndex, SOP>::const_iterator iter = map.find(i);
  assert(iter != map.end());
  return iter->second;
}

}

}
