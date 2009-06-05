#ifndef SIMPLE_CFD_MESH_ENTITY_HPP
#define SIMPLE_CFD_MESH_ENTITY_HPP

namespace cfd
{

class MeshEntity
{
private:
  MeshTopology* topology;
  std::size_t dimension;
  std::size_t index;

public:
  MeshEntity(MeshTopology* const _topology, const std::size_t _dimension, const std::size_t _index);
};

}

#endif
