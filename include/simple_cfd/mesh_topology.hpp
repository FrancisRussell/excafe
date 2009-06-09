#ifndef SIMPLE_CFD_MESH_TOPOLOGY_HPP
#define SIMPLE_CFD_MESH_TOPOLOGY_HPP

#include <mesh_connectivity.hpp>
#include <mesh_entity_iterator_global.hpp>
#include <mesh_entity_iterator_local.hpp>
#include <general_cell.hpp>
#include <cstddef>
#include <cassert>
#include <vector>
#include <set>

namespace cfd
{

class MeshTopology
{
private:
  friend class MeshEntityIteratorGlobal;
  friend class MeshEntityIteratorLocal;

  const GeneralCell& cell;
  const std::size_t dimension;
  std::vector<MeshConnectivity> relations;

  static std::size_t numRelations(const std::size_t dimension);
  std::size_t getConnectivityIndex(const std::size_t d, const std::size_t dPrime) const;
  MeshConnectivity* getConnectivity(const std::size_t d, const std::size_t dPrime);

  void calculateConnectivity(const std::size_t d, const std::size_t dPrime);
  void performTranspose(const std::size_t d, const std::size_t dPrime);
  void performIntersection(const std::size_t d, const std::size_t dPrime, const std::size_t dPrimePrime);
  void performBuild(const std::size_t d);


public:
  typedef MeshEntityIteratorGlobal global_iterator;
  typedef MeshEntityIteratorLocal local_iterator;

  MeshTopology(const GeneralCell& _cell);
  std::size_t numEntities(const std::size_t d);
  std::set<std::size_t> getIndices(const MeshEntity& entity, const std::size_t d);

  global_iterator global_begin(const std::size_t d);
  global_iterator global_end(const std::size_t d);

  local_iterator local_begin(const MeshEntity& entity, const std::size_t d);
  local_iterator local_end(const MeshEntity& entity, const std::size_t d);
};

}

#endif 
