#ifndef SIMPLE_CFD_LAGRANGE_TRIANGLE_LINEAR_HPP
#define SIMPLE_CFD_LAGRANGE_TRIANGLE_LINEAR_HPP

#include <map>
#include <vector>
#include <utility>
#include "simple_cfd_fwd.hpp"
#include "mesh.hpp"
#include "finite_element.hpp"

namespace cfd
{

template<unsigned R>
class lagrange_triangle_linear : public finite_element< cell<triangle> >
{
public:
  typedef cell<triangle> cell_type;
  static const unsigned int rank = R;
  static const unsigned int dimension = cell_type::dimension;
  typedef vertex<dimension> vertex_type;

  const mesh<cell_type>* m;

  lagrange_triangle_linear(const mesh<cell_type>& _m) : m(&_m)
  {
  }

  evaluated_basis evaluate_basis(const cell_type& c, const unsigned int i, const vertex_type& v) const
  {
    const mesh_geometry<dimension> geometry(m->getGeometry());
    const std::vector<vertex_type> vertices(c.getCoordinates(geometry));
    const double area = c.getArea(geometry);

    const int ip1 = (i+1) % 3;
    const int ip2 = (i+2) % 3;

    evaluated_basis result;

    result.value = ((vertices[ip2][0] - vertices[ip1][0]) * (v[1] - vertices[ip1][1]) -
                    (vertices[ip2][1] - vertices[ip1][1]) * (v[0] - vertices[ip1][0])) / (2.0 * area);
    result.dx = -(vertices[ip2][1] - vertices[ip1][1]) / (2.0 * area);
    result.dy =  (vertices[ip2][0] - vertices[ip1][0]) / (2.0 * area);

    return result;
  }

  unsigned space_dimension() const
  {
    return 3;
  }

  std::vector< std::pair<unsigned, unsigned> > getCommonDegreesOfFreedom(const cell_id cid, const cell_id cid2) const
  {
    const std::vector<vertex_id> cid_vertices(m->getCell(cid).getIndices());
    const std::vector<vertex_id> cid2_vertices(m->getCell(cid2).getIndices());

    // Map vertices of cid2 onto degrees of freedom
    std::map<vertex_id, unsigned> cid2_dof;
    for(unsigned dof=0; dof<cid2_vertices.size(); ++dof)
      cid2_dof[cid2_vertices[dof]] = dof;

    std::vector< std::pair<unsigned, unsigned> > common;
    // Iterate over degrees of freedom on cid and find ones that correspond to common vertices
    for(unsigned dof=0; dof<cid_vertices.size(); ++dof)
    {
      const std::map<vertex_id, unsigned>::const_iterator sharedVertexIter = cid2_dof.find(cid_vertices[dof]);

      if (sharedVertexIter != cid2_dof.end())
        common.push_back(std::make_pair(dof, sharedVertexIter->second));
    }
    return common;
  }

  std::vector<unsigned> getBoundaryDegreesOfFreedom(const cell_id cid, const std::vector< std::pair<vertex_id, vertex_id> >& boundary) const
  {
    // Create set of vertices on boundary
    std::set<vertex_id> boundaryVertices;
    for(std::vector< std::pair<vertex_id, vertex_id> >::const_iterator edgeIter(boundary.begin()); edgeIter!=boundary.end(); ++edgeIter)
    {
      boundaryVertices.insert(edgeIter->first);
      boundaryVertices.insert(edgeIter->second);
    }

    // Create list of local vertices
    const std::vector<vertex_id> vertexIndices(m->getCell(cid).getIndices());

    // Find vertices on boundary
    std::vector<unsigned> dofs;
    for(unsigned i=0; i<vertexIndices.size(); ++i)
    {
      if (boundaryVertices.find(vertexIndices[i]) != boundaryVertices.end())
        dofs.push_back(i);
    }
    return dofs;
  }

  vertex_type getDofCoordinate(const cell_id cid, const unsigned dof) const
  {
    assert(dof>=0 && dof<3);
    return m->getCoordinates(cid)[dof];
  }
};

}

#endif
