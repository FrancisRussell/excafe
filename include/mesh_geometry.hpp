#ifndef SIMPLE_CFD_MESH_GEOMETRY_HPP
#define SIMPLE_CFD_MESH_GEOMETRY_HPP

#include "simple_cfd_fwd.hpp"
#include <utility>
#include <cassert>
#include <cstddef>
#include <map>

namespace cfd
{

template<unsigned int D>
class mesh_geometry
{
public:
  static const unsigned int dimension = D;

private:
  std::map<vertex_id, vertex<dimension> > values;

public:
  mesh_geometry()
  {
  }

  std::size_t size() const
  {
    return values.size();
  }

  void insert(const vertex_id vid, const vertex<dimension>& v)
  {
    values.insert(std::make_pair(vid, v));
  }

  vertex<dimension>& operator[](const vertex_id vid)
  {
    return values[vid];
  }

  const vertex<dimension> operator[](const vertex_id vid) const
  {
    const typename std::map<vertex_id, vertex<dimension> >::const_iterator vertexIter(values.find(vid));

    if (vertexIter != values.end())
    {
      return vertexIter->second;
    }
    else
    {
      assert(false);
      return vertex<dimension>();
    }
  }
};

}

#endif
