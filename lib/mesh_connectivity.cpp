#include <mesh_connectivity.hpp>
#include <cstddef>
#include <cassert>
#include <algorithm>

namespace cfd
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

std::size_t MeshConnectivity::entitySize(const std::size_t entity) const
{
  assert(entity < numEntities());
  return offsets[entity+1] - offsets[entity];
}

void MeshConnectivity::populateWithIndices(std::vector<std::size_t>& _indices, const std::size_t entity) const
{
  assert(entity < numEntities());
  _indices.resize(entitySize(entity));
  std::copy(&indices[offsets[entity]], &indices[offsets[entity+1]], _indices.begin());
}

}
