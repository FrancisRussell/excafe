#ifndef SIMPLE_CFD_MESH_TOPOLOGY_HPP
#define SIMPLE_CFD_MESH_TOPOLOGY_HPP

#include <mesh_connectivity.hpp>
#include <mesh_entity_iterator_global.hpp>
#include <cstddef>
#include <cassert>
#include <vector>

namespace cfd
{

class MeshTopology
{
private:
  const std::size_t dimension;
  std::vector<MeshConnectivity> relations;

  static std::size_t numRelations(const std::size_t dimension);
  std::size_t getConnectivityIndex(const std::size_t d, const std::size_t dPrime) const;

  void calculateConnectivity(const std::size_t d, const std::size_t dPrime);
  void performTranspose(const std::size_t d, const std::size_t dPrime);
  void performIntersection(const std::size_t d, const std::size_t dPrime, const std::size_t dPrimePrime);
  void performBuild(const std::size_t d);

public:
  typedef MeshEntityIteratorGlobal global_iterator;

  MeshTopology(const std::size_t _dimension);
  std::size_t numEntities(const std::size_t d);

//  iterator begin(d, dPrime);
//  iterator end(d, dPrime);

//  const_iterator begin(d, dPrime);
//  const_iterator end(d, dPrime);
};

}

#endif 
