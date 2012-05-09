#ifndef EXCAFE_MESH_GEOMETRY_IMPL_HPP
#define EXCAFE_MESH_GEOMETRY_IMPL_HPP

#include "excafe_fwd.hpp"
#include <vector>
#include <cassert>
#include <cstddef>

namespace excafe
{

template<unsigned int D>
class MeshGeometryImpl
{
public:
  static const unsigned int dimension = D;

private:
  std::vector< vertex<dimension> > values;

public:
  MeshGeometryImpl()
  {
  }

  std::size_t size() const
  {
    return values.size();
  }

  const vertex_id add(const vertex<dimension>& v)
  {
    const vertex_id vid = values.size();
    values.push_back(v);
    return vid;
  }

  vertex<dimension>& operator[](const vertex_id vid)
  {
    assert(vid < values.size());
    return values[vid];
  }

  const vertex<dimension> operator[](const vertex_id vid) const
  {
    assert(vid < values.size());
    return values[vid];
  }
};

}

#endif
