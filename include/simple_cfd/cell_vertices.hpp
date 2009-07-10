#ifndef SIMPLE_CFD_CELL_VERTICES_HPP
#define SIMPLE_CFD_CELL_VERTICES_HPP

#include <cstddef>
#include <vector>
#include <cassert>
#include <vertex.hpp>

namespace cfd
{

template<std::size_t D>
class CellVertices
{
public:
  static const std::size_t dimension = D;
  typedef typename std::vector< vertex<dimension> >::const_iterator const_iterator;
  typedef const_iterator iterator;

private:
  std::vector< vertex<dimension> > vertices;

public:
  template<typename InputIterator>
  CellVertices(const InputIterator& begin, const InputIterator& end) : vertices(begin, end)
  {
  }

  const_iterator begin() const
  {
    return vertices.begin();
  }

  const_iterator end() const
  {
    return vertices.end();
  }

  const vertex<dimension> operator[](const std::size_t index) const
  {
    assert(index < vertices.size());
    return vertices[index];
  }
};

}

#endif
