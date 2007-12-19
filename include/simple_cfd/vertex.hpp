#ifndef SIMPLE_CFD_VERTEX_HPP
#define SIMPLE_CFD_VERTEX_HPP

#include "simple_cfd_fwd.hpp"
#include <cassert>
#include <algorithm>
#include <vector>
#include <boost/lambda/lambda.hpp>
#include <boost/array.hpp>
#include <boost/operators.hpp>

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
};

}

#endif
