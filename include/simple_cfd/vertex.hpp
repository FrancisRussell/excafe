#ifndef SIMPLE_CFD_VERTEX_HPP
#define SIMPLE_CFD_VERTEX_HPP

#include "simple_cfd_fwd.hpp"
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include <boost/lambda/lambda.hpp>
#include <boost/array.hpp>
#include <boost/operators.hpp>
#include <ostream>

namespace cfd
{

template<std::size_t D, typename T>
class vertex :  boost::totally_ordered<vertex<D, T>,
                boost::additive<vertex<D, T>,
                boost::multiplicative<vertex<D, T>, T
                > > >
{
public:
  static const std::size_t dimension = D;
  typedef T value_type;
  typedef typename boost::array<value_type, dimension>::iterator iterator;
  typedef typename boost::array<value_type, dimension>::const_iterator const_iterator;

private:
  boost::array<value_type, dimension> values;

public:
  vertex()
  {
  }

  vertex(const vertex& v) : values(v.values)
  {
  }

  vertex(const value_type& x)
  {
    assert(dimension == 1);
    values[0] = x;
  }

  vertex(const value_type& x, const value_type& y)
  {
    assert(dimension == 2);
    values[0] = x;
    values[1] = y;
  }

  vertex(const std::vector<value_type>& v)
  {
    assert(values.size() == dimension);
    std::copy(v.begin(), v.end(), values.begin());
  }

  value_type& operator[](const unsigned int index)
  {
    return values[index];
  }
 
  const value_type& operator[](const unsigned int index) const
  {
    return values[index];
  }
 
  bool operator==(const vertex& v) const
  {
    return std::equal(values.begin(), values.end(), v.values.begin());
  }

  bool operator<(const vertex& v) const
  {
    return std::lexicographical_compare(values.begin(), values.end(), v.values.begin(), v.values.end());
  }

  vertex& operator+=(const vertex& v) 
  {
    std::transform(values.begin(), values.end(), v.values.begin(), values.begin(),
      boost::lambda::ret<value_type>(boost::lambda::_1 + boost::lambda::_2));
    return *this;
  }

  vertex& operator-=(const vertex& v)
  {
    std::transform(values.begin(), values.end(), v.values.begin(), values.begin(), 
      boost::lambda::ret<value_type>(boost::lambda::_1 - boost::lambda::_2));
    return *this;
  }

  vertex& operator/=(const value_type& d)
  {
    std::transform(values.begin(), values.end(), values.begin(), 
      boost::lambda::ret<value_type>(boost::lambda::_1 / d));
    return *this;
  }

  vertex& operator*=(const value_type& d)
  {
    std::transform(values.begin(), values.end(), values.begin(), 
      boost::lambda::ret<value_type>(boost::lambda::_1 * d));
    return *this;
  }

  iterator begin()
  {
    return values.begin();
  }

  const_iterator begin() const
  {
    return values.begin();
  }

  iterator end()
  {
    return values.end();
  }

  const_iterator end() const
  {
    return values.end();
  }

  value_type distance(const vertex<dimension>& v) const 
  {
    const vertex delta = v - *this;
    return std::sqrt(std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0));
  }
};

template<std::size_t D, typename T>
std::ostream& operator<<(std::ostream& o, const vertex<D, T>& v)
{
  o << "(";

  for(typename vertex<D, T>::const_iterator vIter(v.begin()); vIter!=v.end(); ++vIter)
  {
    o << *vIter;

    if (vIter+1 != v.end())
      o << ", ";
  }

  o << ")";
  return o;
}

}

#endif
