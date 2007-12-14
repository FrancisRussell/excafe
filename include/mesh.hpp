#ifndef SIMPLE_CFD_MESH_HPP
#define SIMPLE_CFD_MESH_HPP

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <boost/lambda/lambda.hpp>
#include "simple_cfd_fwd.hpp"
#include "utility.hpp"
#include "cell.hpp"
#include "dof_map.hpp"
#include "matrix.hpp"

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
  int x_size;
  int y_size;
  std::map<vertex_id, vertex_type> vertices;
  std::map<cell_id, cell_type> cells;
  std::map<vertex_id, std::set<cell_id> > vertex_to_cells;

public:
  mesh(const int width, const int height) : x_size(width), y_size(height)
  {
    assert(x_size > 1);
    assert(y_size > 1);

    // Create vertices
    vertex_id vid=0;

    const int x_nodes = x_size;
    const int y_nodes = y_size;

    for(int y=0; y < y_nodes; ++y)
    {
      for(int x=0; x < x_nodes; ++x)
      {
        addVertex(vid, vertex_type(static_cast<double>(x) / (x_nodes-1), static_cast<double>(y) / (y_nodes-1)));
        ++vid;
      }
    }

    /*   Triangle Node Numbering
  
         LL        UR

         2           1     0
                      _____
         |\           \    |
         | \           \   |
         |  \           \  |  
         |   \           \ |
         |    \           \|
         ------           
         0     1           2
    */

    // Now create cells assuming first vertex has id 0
    std::vector<vertex_id> lower_left_vertices;
    lower_left_vertices.push_back(0);
    lower_left_vertices.push_back(1);
    lower_left_vertices.push_back(x_nodes);

    std::vector<vertex_id> upper_right_vertices;
    upper_right_vertices.push_back(1 + x_nodes);
    upper_right_vertices.push_back(x_nodes);
    upper_right_vertices.push_back(1);

    cell_id cid = 0;

    for(int y = 0; y < y_size - 1; ++y)
    {
      for(int x =0; x < x_size - 1; ++x)
      {
        std::vector<vertex_id> offset_lower_left_vertices(lower_left_vertices);
        std::vector<vertex_id> offset_upper_right_vertices(upper_right_vertices);

        const int offset = x_nodes * y + x;

        std::transform(offset_lower_left_vertices.begin(), 
          offset_lower_left_vertices.end(), offset_lower_left_vertices.begin(), boost::lambda::_1 + offset);

        std::transform(offset_upper_right_vertices.begin(), 
          offset_upper_right_vertices.end(), offset_upper_right_vertices.begin(), boost::lambda::_1 + offset);

        addCell(cid, cell_type(offset_lower_left_vertices));
        ++cid;
        addCell(cid, cell_type(offset_upper_right_vertices));
        ++cid;
      }
    }
  }
  
  void addVertex(const vertex_id vid, const vertex_type& v)
  {
    vertices.insert(std::make_pair(vid, v));
  }

  void addCell(const cell_id cid, const cell_type& c)
  {
    cells.insert(std::make_pair(cid, c));
    const std::vector< vertex_id > cell_vertices(c.getIndices());
    for(std::vector<vertex_id>::const_iterator vertexIter(cell_vertices.begin()); vertexIter != cell_vertices.end(); ++vertexIter)
    {
      vertex_to_cells[*vertexIter].insert(cid);
    }
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

  std::set< std::pair<vertex_id, vertex_id> > getEdgeFacets() const
  {
    std::map< std::pair<vertex_id, vertex_id>, unsigned, unordered_pair_compare<vertex_id> > facetCount;
    for(std::map<cell_id, cell_type>::const_iterator cellIter(cells.begin()); cellIter!=cells.end(); ++cellIter)
    {
      const std::set< std::pair<vertex_id, vertex_id> > facets(cellIter->second.getFacets());
      for(std::set< std::pair<vertex_id, vertex_id> >::const_iterator facetIter(facets.begin()); facetIter!=facets.end(); ++facetIter)
        ++facetCount[*facetIter];
    }

    std::set< std::pair<vertex_id, vertex_id> > facets;
    for(std::map< std::pair<vertex_id, vertex_id>, unsigned >::const_iterator facetCountIter(facetCount.begin()); facetCountIter!=facetCount.end(); ++facetCountIter)
    {
      // We only want facets that are not adjacent to other facets
      if (facetCountIter->second == 1)
        facets.insert(facetCountIter->first);
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
    out << "Nodes: " << vertices.size() << std::endl;
    out << "Cells: " << cells.size() << std::endl;
    out << std::endl;

    for(std::map<cell_id, cell_type>::const_iterator cellIter = cells.begin(); cellIter != cells.end(); ++cellIter)
    {
      out << "Cell: " << cellIter->first << std::endl;
      cellIter->second.print(out);
      out << std::endl;
      //printArea(cellIter->second);
    }
  }
  
  cell_type getCell(const cell_id cid) const
  {
    const std::map<cell_id, cell_type>::const_iterator cellIter = cells.find(cid);
    assert(cellIter != cells.end());
    return cellIter->second;
  }

  vertex_type getVertex(const vertex_id vid) const
  {
    const std::map<vertex_id, vertex_type>::const_iterator vertexIter = vertices.find(vid);
    assert(vertexIter != vertices.end());
    return vertexIter->second;
  }

  std::map<cell_id, cell_type> getCells() const
  {
    return cells;
  }

  std::map<vertex_id, vertex_type> getVertices() const
  {
    return vertices;
  }

  virtual ~mesh()
  {
  }
};

}

#endif
