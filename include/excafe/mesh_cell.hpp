#ifndef EXCAFE_MESH_CELL_HPP
#define EXCAFE_MESH_CELL_HPP

#include <set>
#include <cstddef>
#include <memory>
#include "excafe_fwd.hpp"

namespace excafe
{

class MeshCell
{
public:
  virtual std::size_t getDimension() const = 0;
  virtual std::size_t numEntities(std::size_t dimension) const = 0;
  virtual std::set<std::size_t> getIncidentVertices(const MeshEntity& localEntity) const = 0;
  virtual std::set<std::size_t> getIncidentVertices(MeshTopology& topology, const std::size_t cid, const MeshEntity& localEntity) const = 0;
  virtual ~MeshCell() {}
};

}
#endif
