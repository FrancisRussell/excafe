#include <mesh_entity_iterator_local.hpp>

namespace cfd
{

MeshEntityIteratorLocal::MeshEntityIteratorLocal(MeshTopology* const _topology,
  const MeshEntity& _from, const std::size_t _dimensionTo, const std::size_t _offset) : topology(_topology),
  from(_from), dimensionTo(_dimensionTo), offset(_offset), 
  connectivity(topology->getConnectivity(from.getDimension(), dimensionTo))
{
}

}
