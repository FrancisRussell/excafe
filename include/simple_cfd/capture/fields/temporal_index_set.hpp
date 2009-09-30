#ifndef SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_SET_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_TEMPORAL_INDEX_SET_HPP

#include "temporal_index_value.hpp"
#include <set>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <boost/operators.hpp>
#include <boost/iterator/indirect_iterator.hpp>

namespace cfd
{

namespace detail
{

class TemporalIndexSet : boost::addable<TemporalIndexSet, TemporalIndexSet,
                         boost::subtractable<TemporalIndexSet, TemporalIndexSet,
                         boost::addable<TemporalIndexSet, TemporalIndexValue*,
                         boost::subtractable<TemporalIndexSet, TemporalIndexValue*
                         > > > >
{
private:
  typedef std::set<TemporalIndexValue*> index_set_t;
  index_set_t indices;

public:
  typedef boost::indirect_iterator<index_set_t::iterator> iterator;
  typedef boost::indirect_iterator<index_set_t::const_iterator> const_iterator;

  std::size_t size() const
  {
    return indices.size();
  }

  iterator begin()
  {
    return indices.begin();
  }

  iterator end()
  {
    return indices.end();
  }

  const_iterator begin() const
  {
    return indices.begin();
  }

  const_iterator end() const
  {
    return indices.end();
  }

  bool operator==(const TemporalIndexSet& s) const
  {
    return indices == s.indices;
  }
  
  bool operator<(const TemporalIndexSet& s) const
  {
    return indices < s.indices;
  }

  TemporalIndexSet& operator=(const TemporalIndexSet& s)
  {
    indices = s.indices;
    return *this;
  }

  void swap(TemporalIndexSet& s)
  {
    indices.swap(s.indices);
  }

  TemporalIndexSet& operator+=(const TemporalIndexSet& s)
  {
    indices.insert(s.indices.begin(), s.indices.end());
    return *this;
  }

  TemporalIndexSet& operator+=(TemporalIndexValue* const v)
  {
    indices.insert(v);
    return *this;
  }

  TemporalIndexSet& operator-=(const TemporalIndexSet& s)
  {
    index_set_t subtracted;
    std::set_difference(indices.begin(), indices.end(), s.indices.begin(), s.indices.end(),
      std::inserter(subtracted, subtracted.begin()));
    subtracted.swap(indices);
    return *this;
  }

  TemporalIndexSet& operator-=(TemporalIndexValue* const v)
  {
    indices.erase(v);
    return *this;
  }

  bool contains(const TemporalIndexSet& s) const
  {
    return std::includes(indices.begin(), indices.end(), s.indices.begin(), s.indices.end());
  }

  TemporalIndexSet intersection(const TemporalIndexSet& s) const
  {
    TemporalIndexSet result;
    std::set_intersection(indices.begin(), indices.end(), s.indices.begin(), s.indices.end(),
      std::inserter(result.indices, result.indices.begin()));
    return result;
  }
};

}

}

#endif
