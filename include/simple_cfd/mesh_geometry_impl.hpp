#ifndef SIMPLE_CFD_MESH_GEOMETRY_IMPL_HPP
#define SIMPLE_CFD_MESH_GEOMETRY_IMPL_HPP

#include "simple_cfd_fwd.hpp"
#include <vector>
#include <cassert>
#include <cstddef>

namespace cfd
{

template<unsigned int D>
class mesh_geometry_impl
{
public:
  static const unsigned int dimension = D;

private:
  std::vector< vertex<dimension> > values;

public:
  mesh_geometry_impl()
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
