#ifndef SIMPLE_CFD_MESH_HPP
#define SIMPLE_CFD_MESH_HPP

#include <map>
#include <set>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>
#include "simple_cfd_fwd.hpp"
#include "mesh_geometry.hpp"
#include "mesh_connectivity.hpp"
#include "mesh_topology.hpp"
#include "utility.hpp"
#include "cell.hpp"
#include "dof_map.hpp"

namespace cfd
{

template<typename C>
class mesh
{
};

template<>
class mesh<TriangularCell>
{
public:
  typedef TriangularCell cell_type;
  static const unsigned int dimension = cell_type::dimension;
  typedef vertex<dimension> vertex_type;
  typedef MeshTopology::global_iterator global_iterator;
  typedef MeshTopology::local_iterator local_iterator;

private:
  double width;
  double height;
  mesh_geometry<dimension> geometry;
  MeshConnectivity baseConnectivity;
  TriangularCell referenceCell;
  mutable MeshTopology topology;

public:
  mesh() : topology(referenceCell)
  {
  }

  std::size_t getDimension() const
  {
    return referenceCell.getDimension();
  }

  const vertex_id addVertex(const vertex_type& v)
  {
    return geometry.add(v);
  }

  const cell_id addCell(const std::vector<std::size_t> vertexIndices)
  {
    assert(vertexIndices.size() == referenceCell.getVerticesPerCell());
    const cell_id cid = baseConnectivity.addEntity(vertexIndices.begin(), vertexIndices.end());
    return cid;
  }

  const GeneralCell& getReferenceCell() const
  {
    return referenceCell;
  }

  std::map<vertex_type, double> getQuadrature(const MeshEntity& entity) const
  {
    return referenceCell.getQuadrature(*this, entity);
  }

  double getArea(const std::size_t cid) const
  {
    return referenceCell.getArea(*this, MeshEntity(dimension, cid));
  }

  void finish()
  {
    topology.setBaseConnectivity(baseConnectivity);
    baseConnectivity.clear();
  }

  global_iterator global_begin(const std::size_t d) const
  {
    return topology.global_begin(d);
  }

  global_iterator global_end(const std::size_t d) const
  {
    return topology.global_end(d);
  }

  local_iterator local_begin(const MeshEntity& entity, const std::size_t d) const
  {
    return topology.local_begin(entity, d);
  }

  local_iterator local_end(const MeshEntity& entity, const std::size_t d) const
  {
    return topology.local_end(entity, d);
  }

  std::set<cell_id> getCellIncidentCells(const cell_id cid) const
  {
    const MeshEntity cellEntity(dimension, cid);
    const std::vector<std::size_t> vertices(topology.getIndices(cellEntity, dimension));
    return std::set<std::size_t>(vertices.begin(), vertices.end());
  }

  //NOTE: assumes a 2D mesh
  std::vector< std::pair<vertex_id, vertex_id> > getEdgeFacets() const
  {
    std::vector< std::pair<vertex_id, vertex_id> > result;

    for(MeshTopology::global_iterator facetIter(topology.global_begin(dimension-1)); facetIter!=topology.global_end(dimension-1); ++facetIter)
    {
      if (topology.numRelations(*facetIter, dimension) == 1)
      {
        const std::vector<vertex_id> vertices(topology.getIndices(*facetIter, 0));
        assert(vertices.size() == 2);
        result.push_back(std::make_pair(vertices[0], vertices[1]));
      }
    }
    return result;
  }

  std::vector< vertex<dimension> > getCoordinates(const cell_id cid) const
  {
    const std::vector<vertex_id> vertex_ids(topology.getIndices(MeshEntity(dimension, cid), 0));
    std::vector< vertex<dimension> > coords;

    for(std::vector<vertex_id>::const_iterator vertexIter(vertex_ids.begin()); vertexIter!=vertex_ids.end(); ++vertexIter)
      coords.push_back(getVertex(*vertexIter));

    return coords;
  }

  std::vector<std::size_t> getIndices(const MeshEntity& entity, const std::size_t d) const
  {
    return topology.getIndices(entity, d);
  }

  vertex_type getVertex(const vertex_id vid) const
  {
    return geometry[vid];
  }

  mesh_geometry<dimension> getGeometry() const
  {
    return geometry;
  }

  virtual ~mesh()
  {
  }
};

}

#endif
