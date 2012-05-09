#include <mesh_connectivity.hpp>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <vector>

namespace excafe
{

MeshConnectivity::MeshConnectivity()
{
  clear();
}

void MeshConnectivity::clear()
{
  offsets.clear();
  indices.clear();
  offsets.push_back(0);
}

std::size_t MeshConnectivity::numEntities() const
{
  return offsets.size() - 1;
}

std::size_t MeshConnectivity::addEntity(const std::vector<std::size_t>& _indices)
{
  return addEntity(_indices.begin(), _indices.end());
}

std::size_t MeshConnectivity::numRelations(const std::size_t entity) const
{
  assert(entity < numEntities());
  return offsets[entity+1] - offsets[entity];
}

std::size_t MeshConnectivity::numRelations() const
{
  return indices.size();
}

std::vector<std::size_t> MeshConnectivity::getIndices(const std::size_t entity) const
{
  assert(entity < numEntities());
  const std::vector<std::size_t> result(&indices[offsets[entity]], &indices[offsets[entity+1]]);
  return result;
}

}
