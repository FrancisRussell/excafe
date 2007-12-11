#ifndef SIMPLE_CFD_VERTEX_HPP
#define SIMPLE_CFD_VERTEX_HPP

#include "simple_cfd_fwd.hpp"
#include <cassert>
#include <algorithm>
#include <vector>
#include <boost/lambda/lambda.hpp>
#include <boost/array.hpp>
#include <boost/static_assert.hpp>

namespace cfd
{

template<unsigned int D>
class vertex
{
public:
  static const unsigned int dimension = D;

private:
  boost::array<double, D> values;

public:
  vertex()
  {
  }

  vertex(const double x, const double y)
  {
    BOOST_STATIC_ASSERT(dimension == 2);
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

  vertex operator+(const vertex& v) const
  {
    std::vector<double> new_values(values.begin(), values.end());
    std::transform(new_values.begin(), new_values.end(), v.values.begin(), new_values.begin(), boost::lambda::_1 + boost::lambda::_2);
    return vertex(new_values);
  }

  vertex operator/(const double d) const
  {
    std::vector<double> new_values(values.begin(), values.end());
    std::transform(new_values.begin(), new_values.end(), new_values.begin(), boost::lambda::_1 / d);
    return vertex(new_values);
  }

};

}

#endif
