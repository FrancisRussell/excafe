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
class mesh_geometry
{
public:
  static const unsigned int dimension = D;

private:
  boost::shared_ptr< mesh_geometry_impl<D> > impl;

  void make_unique()
  {
    if (!impl.unique())
      impl =  boost::shared_ptr< mesh_geometry_impl<D> >(new mesh_geometry_impl<D>(*impl));
  }

public:
  mesh_geometry() : impl(new mesh_geometry_impl<D>())
  {
  }

  std::size_t size() const
  {
    return impl->size();
  }

  void insert(const vertex_id vid, const vertex<dimension>& v)
  {
    make_unique();
    impl->insert(vid, v);
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
