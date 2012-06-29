#ifndef EXCAFE_CSE_VERTEX_INFO_HPP
#define EXCAFE_CSE_VERTEX_INFO_HPP

#include <boost/operators.hpp>
#include <cstddef>

namespace excafe
{

namespace cse
{

template<typename G>
class VertexInfo : boost::addable< VertexInfo<G> >
{
private:
  std::size_t count;
  std::size_t numeric;
  std::size_t unit;
  std::size_t haveCoefficients;
  int value;

public:
  typedef G graph_t;
  typedef typename graph_t::vertex_descriptor vertex_descriptor;

  VertexInfo() : count(0), numeric(0), unit(0), haveCoefficients(0), value(0)
  {
  }

  void addVertex(const graph_t& graph, const vertex_descriptor& v)
  {
    ++count;
    const typename boost::vertex_bundle_type<graph_t>::type& vertexProperty = graph[v];
    numeric += (vertexProperty.is_numeric ? 1 : 0);
    unit += (vertexProperty.is_unit ? 1 : 0);
    value += vertexProperty.mul_count;
    haveCoefficients += (vertexProperty.has_coefficient ? 1 : 0);
  }

  VertexInfo& operator+=(const VertexInfo& v)
  {
    count += v.count;
    numeric += v.numeric;
    unit += v.unit;
    value += v.value;
    return *this;
  }

  std::size_t num() const
  {
    return count;
  }

  std::size_t numNumeric() const
  {
    return numeric;
  }

  std::size_t numUnit() const
  {
    return unit;
  }

  std::size_t numNonUnit() const
  {
    return count - unit;
  }

  std::size_t numHaveCoefficients() const
  {
    return haveCoefficients;
  }

  int getValue() const
  {
    return value;
  }
};

}

}

#endif
