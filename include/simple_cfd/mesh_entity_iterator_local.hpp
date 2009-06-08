#ifndef SIMPLE_CFD_MESH_ENTITY_ITERATOR_LOCAL
#define SIMPLE_CFD_MESH_ENTITY_ITERATOR_LOCAL

#include <cstddef>
#include <boost/iterator/iterator_facade.hpp>

namespace cfd
{

class MeshEntityIteratorLocal : public boost::iterator_facade<MeshEntityIteratorLocal, const MeshEntity, boost::random_access_iterator_tag>
{
private:
  MeshTopology* topology;
  std::size_t dimension;
  std::size_t index;

public:
  const MeshEntity& dereference() const
  {
  }

  bool equal(const MeshEntityIteratorLocal& e) const
  {
  }

  void increment()
  {
  }

  void decrement()
  {
  }

  void advance(const std::size_t n)
  {
  }

  std::size_t distance_to(const MeshEntityIteratorLocal& e) const
  {
  }
};


}

#endif
