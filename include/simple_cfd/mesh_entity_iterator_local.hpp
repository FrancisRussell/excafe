#ifndef SIMPLE_CFD_MESH_ENTITY_ITERATOR_LOCAL
#define SIMPLE_CFD_MESH_ENTITY_ITERATOR_LOCAL

#include <cstddef>
#include <iterator>
#include <mesh_topology.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace cfd
{

class MeshEntityIteratorLocal : public boost::iterator_facade<MeshEntityIteratorLocal, const MeshEntity,
  std::random_access_iterator_tag, const MeshEntity>
{
private:
  MeshTopology* topology;
  MeshEntity from;
  std::size_t dimensionTo;
  std::size_t offset;

  MeshConnectivity* connectivity;

public:
  MeshEntityIteratorLocal(MeshTopology* const _topology,
    const MeshEntity& _from, const std::size_t _dimensionTo, const std::size_t _offset);

  const MeshEntity dereference() const
  {
    // FIXME: I'm hideously inefficient
    std::vector<std::size_t> indices(connectivity->getIndices(from.getIndex()));
    return MeshEntity(dimensionTo, indices[offset]);
  }

  bool equal(const MeshEntityIteratorLocal& e) const
  {
    return topology == e.topology &&
      from == e.from &&
      dimensionTo == e.dimensionTo &&
      offset == e.offset;
  }

  void increment()
  {
    ++offset;
  }

  void decrement()
  {
    --offset;
  }

  void advance(const std::size_t n)
  {
    offset += n;
  }

  std::size_t distance_to(const MeshEntityIteratorLocal& e) const
  {
    return e.offset - offset;
  }
};


}

#endif
