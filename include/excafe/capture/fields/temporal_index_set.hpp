#ifndef EXCAFE_CAPTURE_FIELDS_TEMPORAL_INDEX_SET_HPP
#define EXCAFE_CAPTURE_FIELDS_TEMPORAL_INDEX_SET_HPP

#include "fields_fwd.hpp"
#include <set>
#include <cstddef>
#include <boost/operators.hpp>
#include <boost/iterator/indirect_iterator.hpp>

namespace excafe
{

namespace detail
{

class TemporalIndexSet : boost::additive<TemporalIndexSet, TemporalIndexSet,
                         boost::additive<TemporalIndexSet, TemporalIndexValue*
                         > >
{
private:
  typedef std::set<TemporalIndexValue*> index_set_t;
  index_set_t indices;

public:
  typedef boost::indirect_iterator<index_set_t::iterator> iterator;
  typedef boost::indirect_iterator<index_set_t::const_iterator> const_iterator;

  std::size_t size() const;
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;
  bool operator==(const TemporalIndexSet& s) const;
  bool operator<(const TemporalIndexSet& s) const;
  TemporalIndexSet& operator=(const TemporalIndexSet& s);
  void swap(TemporalIndexSet& s);
  TemporalIndexSet& operator+=(const TemporalIndexSet& s);
  TemporalIndexSet& operator+=(TemporalIndexValue* const v);
  TemporalIndexSet& operator-=(const TemporalIndexSet& s);
  TemporalIndexSet& operator-=(TemporalIndexValue* const v);
  bool contains(const TemporalIndexSet& s) const;
  TemporalIndexSet intersection(const TemporalIndexSet& s) const;
};

}

}

#endif
