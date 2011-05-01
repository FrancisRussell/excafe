#ifndef SIMPLE_CFD_NUMERIC_ORTHOTOPE_HPP
#define SIMPLE_CFD_NUMERIC_ORTHOTOPE_HPP

#include <map>
#include <set>
#include <vector>
#include <utility>
#include <boost/foreach.hpp>

namespace cfd
{

namespace numeric
{

template<typename V, typename N>
class Orthotope
{
private:
  typedef V variable_t;
  typedef N numeric_t;

  typedef std::vector<std::pair<variable_t, std::pair<numeric_t, numeric_t> > > range_map_t;

  range_map_t ranges;

public:
  typedef typename range_map_t::iterator       iterator;
  typedef typename range_map_t::const_iterator const_iterator;
  typedef typename range_map_t::value_type     value_type;

  Orthotope()
  {
  }

  iterator begin() 
  {
    return ranges.begin();
  }

  const_iterator begin() const 
  {
    return ranges.begin();
  }

  iterator end()
  {
    return ranges.end();
  }

  const_iterator end() const
  {
    return ranges.end();
  }

  std::set<variable_t> getVariables()
  {
    std::set<variable_t> result;
    BOOST_FOREACH(const typename range_map_t::value_type& range, ranges)
      result.insert(range.first);

    return result;
  }

  void setInterval(const variable_t& v, const numeric_t lower, const numeric_t upper)
  {
    ranges.push_back(std::make_pair(v, std::make_pair(lower, upper)));
  }

  numeric_t getVolume() const
  {
    numeric_t result = 1;
    BOOST_FOREACH(const typename range_map_t::value_type& range, ranges)
  }
};

}

}

#endif
