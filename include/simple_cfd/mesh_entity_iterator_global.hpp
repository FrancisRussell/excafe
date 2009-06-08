#ifndef SIMPLE_CFD_MESH_ENTITY_ITERATOR_GLOBAL
#define SIMPLE_CFD_MESH_ENTITY_ITERATOR_GLOBAL

#include <cstddef>
#include <iterator>
#include <simple_cfd_fwd.hpp>
#include <mesh_entity.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace cfd
{

class MeshEntityIteratorGlobal : public boost::iterator_facade<MeshEntityIteratorGlobal, const MeshEntity,
  std::random_access_iterator_tag, const MeshEntity>
{
private:
  MeshTopology* topology;
  std::size_t dimension;
  std::size_t index;

public:
  MeshEntityIteratorGlobal(MeshTopology* const _topology, const std::size_t _dimension) : topology(_topology),
    dimension(_dimension), index(0)
  {
  }

  const MeshEntity dereference() const
  {
    return MeshEntity(dimension, index);
  }

  bool equal(const MeshEntityIteratorGlobal& e) const
  {
    return topology == e.topology &&
      dimension == e.dimension &&
      index == e.index;
  }

  void increment()
  {
    ++index;
  }

  void decrement()
  {
    --index;
  }

  void advance(const std::size_t n)
  {
    index += n;
  }

  std::size_t distance_to(const MeshEntityIteratorGlobal& e) const
  {
    return e.index - index;
  }
};


}

#endif
