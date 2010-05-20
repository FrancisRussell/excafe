#include <cstddef>
#include <iterator>
#include <algorithm>
#include <simple_cfd/capture/fields/temporal_index_set.hpp>
#include <simple_cfd/capture/fields/temporal_index_value.hpp>

namespace cfd
{

namespace detail
{

std::size_t TemporalIndexSet::size() const
{
  return indices.size();
}

TemporalIndexSet::iterator TemporalIndexSet::begin()
{
  return indices.begin();
}

TemporalIndexSet::iterator TemporalIndexSet::end()
{
  return indices.end();
}

TemporalIndexSet::const_iterator TemporalIndexSet::begin() const
{
  return indices.begin();
}

TemporalIndexSet::const_iterator TemporalIndexSet::end() const
{
  return indices.end();
}

bool TemporalIndexSet::operator==(const TemporalIndexSet& s) const
{
  return indices == s.indices;
}

bool TemporalIndexSet::operator<(const TemporalIndexSet& s) const
{
  return indices < s.indices;
}

TemporalIndexSet& TemporalIndexSet::operator=(const TemporalIndexSet& s)
{
  indices = s.indices;
  return *this;
}

void TemporalIndexSet::swap(TemporalIndexSet& s)
{
  indices.swap(s.indices);
}

TemporalIndexSet& TemporalIndexSet::operator+=(const TemporalIndexSet& s)
{
  indices.insert(s.indices.begin(), s.indices.end());
  return *this;
}

TemporalIndexSet& TemporalIndexSet::operator+=(TemporalIndexValue* const v)
{
  indices.insert(v);
  return *this;
}

TemporalIndexSet& TemporalIndexSet::operator-=(const TemporalIndexSet& s)
{
  index_set_t subtracted;
  std::set_difference(indices.begin(), indices.end(), s.indices.begin(), s.indices.end(),
    std::inserter(subtracted, subtracted.begin()));
  subtracted.swap(indices);
  return *this;
}

TemporalIndexSet& TemporalIndexSet::operator-=(TemporalIndexValue* const v)
{
  indices.erase(v);
  return *this;
}

bool TemporalIndexSet::contains(const TemporalIndexSet& s) const
{
  return std::includes(indices.begin(), indices.end(), s.indices.begin(), s.indices.end());
}

TemporalIndexSet TemporalIndexSet::intersection(const TemporalIndexSet& s) const
{
  TemporalIndexSet result;
  std::set_intersection(indices.begin(), indices.end(), s.indices.begin(), s.indices.end(),
    std::inserter(result.indices, result.indices.begin()));
  return result;
}

}

}
