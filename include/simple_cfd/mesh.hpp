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
class mesh< cell<triangle> >
{
public:
  static const shape cell_shape = triangle;
  typedef cell<cell_shape> cell_type;
  static const unsigned int dimension = cell_type::dimension;
  typedef vertex<dimension> vertex_type;

private:
  double width;
  double height;
  mesh_geometry<dimension> geometry;
  MeshConnectivity baseConnectivity;
  cell<triangle> referenceCell;
  mutable MeshTopology topology;

public:
  mesh() : topology(referenceCell)
  {
  }

  const vertex_id addVertex(const vertex_type& v)
  {
    return geometry.add(v);
  }

  const cell_id addCell(const cell_type& c)
  {
    const std::vector<vertex_id> cell_vertices(c.getIndices());
    const cell_id cid = baseConnectivity.addEntity(cell_vertices.begin(), cell_vertices.end());
    return cid;
  }

  void finish()
  {
    topology.setBaseConnectivity(baseConnectivity);
    baseConnectivity.clear();
  }

  std::set<cell_id> getCellIncidentCells(const cell_id cid) const
  {
    const MeshEntity cellEntity(dimension, cid);
    return topology.getIndices(cellEntity, dimension);
  }

  //NOTE: assumes a 2D mesh
  std::vector< std::pair<vertex_id, vertex_id> > getEdgeFacets() const
  {
    std::vector< std::pair<vertex_id, vertex_id> > result;

    for(MeshTopology::global_iterator facetIter(topology.global_begin(dimension-1)); facetIter!=topology.global_end(dimension-1); ++facetIter)
    {
      if (topology.numRelations(*facetIter, dimension) == 1)
      {
        const std::set<vertex_id> vertexSet(topology.getIndices(*facetIter, 0));
        const std::vector<vertex_id> vertices(vertexSet.begin(), vertexSet.end());
        assert(vertices.size() == 2);
        result.push_back(std::make_pair(vertices[0], vertices[1]));
      }
    }
    return result;
  }

  std::vector< vertex<dimension> > getCoordinates(const cell_id cid) const
  {
    const std::vector<vertex_id> vertex_ids(getCell(cid).getIndices());
    std::vector< vertex<dimension> > coords;

    for(std::vector<vertex_id>::const_iterator vertexIter(vertex_ids.begin()); vertexIter!=vertex_ids.end(); ++vertexIter)
      coords.push_back(getVertex(*vertexIter));

    return coords;
  }

  void print(std::ostream& out = std::cout) const
  {
    const std::map<cell_id, cell_type> cells(getCells());
    out << "Nodes: " << geometry.size() << std::endl;
    out << "Cells: " << cells.size() << std::endl;
    out << std::endl;

    for(std::map<cell_id, cell_type>::const_iterator cellIter = cells.begin(); cellIter != cells.end(); ++cellIter)
    {
      out << "Cell: " << cellIter->first << std::endl;
      cellIter->second.print(out);
      out << std::endl;
    }
  }
  
  cell_type getCell(const cell_id cid) const
  {
    std::vector<std::size_t> vertexIndices;
    topology.getConnectivity(dimension, 0)->populateWithIndices(vertexIndices, cid);
    return cell_type(vertexIndices);
  }

  vertex_type getVertex(const vertex_id vid) const
  {
    return geometry[vid];
  }

  std::map<cell_id, cell_type> getCells() const
  {
    std::vector<std::size_t> vertexIndices;
    std::map<cell_id, cell_type> cells;

    for(std::size_t cid = 0; cid < topology.numEntities(dimension); ++cid)
    {
      topology.getConnectivity(dimension, 0)->populateWithIndices(vertexIndices, cid);
      cell_type cell(vertexIndices);
      cells.insert(std::make_pair(cid, cell));
    }
    return cells;
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
