#ifndef EXCAFE_MESH_TOPOLOGY_HPP
#define EXCAFE_MESH_TOPOLOGY_HPP

#include <cstddef>
#include <cassert>
#include <vector>
#include <set>
#include <map>
#include <boost/scoped_ptr.hpp>
#include "mesh_connectivity.hpp"
#include "mesh_entity_iterator_global.hpp"
#include "mesh_entity_iterator_local.hpp"
#include "mesh_cell.hpp"
#include "cell_manager.hpp"

namespace excafe
{

class MeshTopology
{
private:
  friend class MeshEntityIteratorGlobal;
  friend class MeshEntityIteratorLocal;

  CellManager::mesh_cell_ref cell;
  const std::size_t dimension;
  std::vector<MeshConnectivity> relations;

  static std::size_t numConnectivityRelations(const std::size_t dimension);
  std::size_t getConnectivityIndex(const std::size_t d, const std::size_t dPrime) const;
  MeshConnectivity* getConnectivityObject(const std::size_t d, const std::size_t dPrime);

  void calculateConnectivity(const std::size_t d, const std::size_t dPrime);
  void performTranspose(const std::size_t d, const std::size_t dPrime);
  void performIntersection(const std::size_t d, const std::size_t dPrime, const std::size_t dPrimePrime);
  void performBuild(const std::size_t d);
  void performBuildZeroToZero();
  bool contains(const MeshEntity& m1, const MeshEntity& m2);

public:
  typedef MeshEntityIteratorGlobal global_iterator;
  typedef MeshEntityIteratorLocal local_iterator;

  MeshTopology(const CellManager::mesh_cell_ref _cell);
  MeshTopology(const MeshTopology& t);

  void setBaseConnectivity(const MeshConnectivity& connectivity);
  std::size_t numEntities(const std::size_t d);
  std::size_t numRelations(const MeshEntity& entity, const std::size_t d);
  std::size_t numRelations(const std::size_t d, const std::size_t dPrime);
  std::vector<std::size_t> getIndices(const MeshEntity& entity, const std::size_t d);

  template<typename OutputIterator>
  void outputIndices(const MeshEntity& entity, const std::size_t d, OutputIterator out)
  {
    const std::vector<std::size_t> indices(getIndices(entity, d));
    std::copy(indices.begin(), indices.end(), out);

  }

  MeshConnectivity* getConnectivity(const std::size_t d, const std::size_t dPrime);

  global_iterator global_begin(const std::size_t d);
  global_iterator global_end(const std::size_t d);

  local_iterator local_begin(const MeshEntity& entity, const std::size_t d);
  local_iterator local_end(const MeshEntity& entity, const std::size_t d);
};

}

#endif 
