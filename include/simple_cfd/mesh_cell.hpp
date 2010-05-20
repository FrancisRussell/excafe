#ifndef SIMPLE_CFD_MESH_CELL_HPP
#define SIMPLE_CFD_MESH_CELL_HPP

#include <vector>
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
  virtual std::vector< std::set<std::size_t> > getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, const std::size_t d) const = 0;
  virtual std::auto_ptr<MeshCell> cloneMeshCell() const = 0;
  virtual ~MeshCell() {}
};

}
#endif
