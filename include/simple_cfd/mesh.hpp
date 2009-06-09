#ifndef SIMPLE_CFD_MESH_HPP
#define SIMPLE_CFD_MESH_HPP

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <boost/lambda/lambda.hpp>
#include "simple_cfd_fwd.hpp"
#include "mesh_geometry.hpp"
#include "mesh_connectivity.hpp"
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
  std::map<vertex_id, std::set<cell_id> > vertex_to_cells;
  cell<triangle> referenceCell;

public:
  mesh()
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
    for(std::vector<vertex_id>::const_iterator vertexIter(cell_vertices.begin()); vertexIter != cell_vertices.end(); ++vertexIter)
    {
      vertex_to_cells[*vertexIter].insert(cid);
    }
    return cid;
  }
  
  std::set<cell_id> getVertexIncidentCells(const vertex_id vid) const
  {
    const std::map<vertex_id, std::set<cell_id> >::const_iterator vertexIter(vertex_to_cells.find(vid));
    assert(vertexIter != vertex_to_cells.end());
    return vertexIter->second;
  }

  std::set<cell_id> getCellIncidentCells(const cell_id cid) const
  {
    // Note that the returned list will also contain the original cell
    const std::vector<vertex_id> indices(getCell(cid).getIndices());
    std::set<cell_id> incident;
    
    for(std::vector<vertex_id>::const_iterator vertexIter(indices.begin()); vertexIter != indices.end(); ++vertexIter)
    {
      const std::set<cell_id> localIncident(getVertexIncidentCells(*vertexIter));
      incident.insert(localIncident.begin(), localIncident.end());
    }
    return incident;
  }

  std::vector< std::pair<vertex_id, vertex_id> > getEdgeFacets() const
  {
    const std::map<cell_id, cell_type> cells(getCells());
    std::map< std::pair<vertex_id, vertex_id>, unsigned, unordered_pair_compare<vertex_id> > facetCount;
    for(std::map<cell_id, cell_type>::const_iterator cellIter(cells.begin()); cellIter!=cells.end(); ++cellIter)
    {
      const std::set< std::pair<vertex_id, vertex_id> > facets(cellIter->second.getFacets());
      for(std::set< std::pair<vertex_id, vertex_id> >::const_iterator facetIter(facets.begin()); facetIter!=facets.end(); ++facetIter)
        ++facetCount[*facetIter];
    }

    std::vector< std::pair<vertex_id, vertex_id> > facets;
    for(std::map< std::pair<vertex_id, vertex_id>, unsigned >::const_iterator facetCountIter(facetCount.begin()); facetCountIter!=facetCount.end(); ++facetCountIter)
    {
      // We only want facets that are not adjacent to other facets
      if (facetCountIter->second == 1)
        facets.push_back(facetCountIter->first);
    }

    return facets;
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
    baseConnectivity.populateWithIndices(vertexIndices, cid);
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

    for(std::size_t cid = 0; cid < baseConnectivity.numEntities(); ++cid)
    {
      baseConnectivity.populateWithIndices(vertexIndices, cid);
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
