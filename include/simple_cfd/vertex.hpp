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

template<unsigned int D>
class vertex :  private boost::less_than_comparable< vertex<D> >,
                private boost::equality_comparable< vertex<D> >,
                private boost::addable< vertex<D> >,
                private boost::subtractable< vertex<D> >,
                private boost::dividable2<vertex<D>, double>,
                private boost::multipliable2<vertex<D>, double>
{
public:
  static const unsigned int dimension = D;
  typedef typename boost::array<double, dimension>::iterator iterator;
  typedef typename boost::array<double, dimension>::const_iterator const_iterator;

private:
  boost::array<double, dimension> values;

public:
  vertex()
  {
  }

  vertex(const vertex& v) : values(v.values)
  {
  }

  vertex(const double x)
  {
    assert(dimension == 1);
    values[0] = x;
  }

  vertex(const double x, const double y)
  {
    assert(dimension == 2);
    values[0] = x;
    values[1] = y;
  }

  vertex(const std::vector<double>& v)
  {
    assert(values.size() == dimension);
    std::copy(v.begin(), v.end(), values.begin());
  }

  double& operator[](const unsigned int index)
  {
    return values[index];
  }
 
  const double& operator[](const unsigned int index) const
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
    std::transform(values.begin(), values.end(), v.values.begin(), values.begin(), boost::lambda::_1 + boost::lambda::_2);
    return *this;
  }

  vertex& operator-=(const vertex& v)
  {
    std::transform(values.begin(), values.end(), v.values.begin(), values.begin(), boost::lambda::_1 - boost::lambda::_2);
    return *this;
  }

  vertex& operator/=(const double d)
  {
    std::transform(values.begin(), values.end(), values.begin(), boost::lambda::_1 / d);
    return *this;
  }

  vertex& operator*=(const double d)
  {
    std::transform(values.begin(), values.end(), values.begin(), boost::lambda::_1 * d);
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

  double distance(const vertex<dimension>& v) const 
  {
    const vertex delta = v - *this;
    return std::sqrt(std::inner_product(delta.begin(), delta.end(), delta.begin(), 0.0));
  }
};

template<unsigned int D>
std::ostream& operator<<(std::ostream& o, const vertex<D>& v)
{
  o << "(";

  for(typename vertex<D>::const_iterator vIter(v.begin()); vIter!=v.end(); ++vIter)
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
