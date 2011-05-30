#ifndef SIMPLE_CFD_MESH_CELL_HPP
#define SIMPLE_CFD_MESH_CELL_HPP

#include <set>
#include <cstddef>
#include <memory>
#include "simple_cfd_fwd.hpp"

namespace cfd
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
