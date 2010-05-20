#ifndef SIMPLE_CFD_MESH_GEOMETRY_HPP
#define SIMPLE_CFD_MESH_GEOMETRY_HPP

#include "simple_cfd_fwd.hpp"
#include <cassert>
#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "mesh_geometry_impl.hpp"

namespace cfd
{

template<unsigned int D>
class MeshGeometry
{
public:
  static const unsigned int dimension = D;

private:
  typedef MeshGeometryImpl<D> impl_t;
  boost::shared_ptr<impl_t> impl;

  void make_unique()
  {
    if (!impl.unique())
      impl =  boost::shared_ptr<impl_t>(new impl_t(*impl));
  }

public:
  MeshGeometry() : impl(new impl_t())
  {
  }

  std::size_t size() const
  {
    return impl->size();
  }

  vertex_id add(const vertex<dimension>& v)
  {
    make_unique();
    return impl->add(v);
  }

  vertex<dimension>& operator[](const vertex_id vid)
  {
    make_unique();
    return (*impl)[vid];
  }

  const vertex<dimension> operator[](const vertex_id vid) const
  {
    return (*impl)[vid];
  }
};

}

#endif
