#ifndef SIMPLE_CFD_LAGRANGE_TRIANGLE_QUADRATIC_HPP
#define SIMPLE_CFD_LAGRANGE_TRIANGLE_QUADRATIC_HPP

#include <vector>
#include <utility>
#include <map>
#include "utility.hpp"
#include "simple_cfd_fwd.hpp"
#include "finite_element.hpp"

namespace cfd
{

class lagrange_triangle_quadratic : public finite_element< cell<triangle> >
{
public:
  typedef cell<triangle> cell_type;
  static const unsigned int dimension = cell_type::dimension;
  typedef vertex<dimension> vertex_type;

  const mesh<cell_type>* m;

  lagrange_triangle_quadratic(const mesh<cell_type>& _m) : m(&_m)
  {
  }

  evaluated_basis evaluate_basis(const cell_id cid, const unsigned int i, const vertex_type& v) const
  {
    /*   Triangle Node Numbering for quadratic basis
  
         LL        

         2
 
         |\ 
         | \ 
       5 |  \ 4
         |   \ 
         |    \ 
         ------ 
         0  3  1
    */

    std::vector<vertex_type> vertices(m->getCoordinates(cid));

    // Create interpolated vertices
    vertices.push_back((vertices[0] + vertices[1])/2);
    vertices.push_back((vertices[1] + vertices[2])/2);
    vertices.push_back((vertices[2] + vertices[0])/2);

    int j1, j2, k1, k2;
    if (i < 3)
    {
      j1 = (i+1)%3;
      j2 = (i+2)%3;
      k1 = 3 + i;
      k2 = 3 + (i + 5)%3;
    }
    else
    {
      j1 = i-3;
      j2 = (i-3+2)%3;
      k1 = (i-3+1)%3;
      k2 = (i-3+2)%3;
    }

    const double gf = (v[0] - vertices[j1][1]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (v[1] - vertices[j1][1]);

    const double gn = (vertices[i][0] - vertices[j1][0]) * (vertices[j2][1] - vertices[j1][1]) -
                      (vertices[j2][0] - vertices[j1][0]) * (vertices[i][1] - vertices[j1][1]);

    const double hf = (v[0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (v[1] - vertices[k1][1]);

    const double hn = (vertices[i][0] - vertices[k1][0]) * (vertices[k2][1] - vertices[k1][1]) -
                      (vertices[k2][0] - vertices[k1][0]) * (vertices[i][1] - vertices[k1][1]);

    evaluated_basis result;

    result.value = (gf/gn) * (hf/hn);

    result.dx = ((vertices[j2][1] - vertices[j1][1])/gn) * (hf/hn) +
                (gf/gn) * ((vertices[k2][1] - vertices[k1][1])/hn);

    result.dy = -((vertices[j2][0] - vertices[j1][0])/gn) * (hf/hn) -
                 (gf/gn) * ((vertices[k2][0] - vertices[k1][0])/hn);

    return result;
  }

  unsigned space_dimension() const
  {
    return 6;
  }

  std::vector< std::pair<unsigned, unsigned> > getCommonDegreesOfFreedom(const cell_id cid, const cell_id cid2) const
  {
    const std::vector<vertex_id> cid_vertices(m->getCell(cid).getIndices());
    const std::vector<vertex_id> cid2_vertices(m->getCell(cid2).getIndices());

    std::vector< std::pair<vertex_id, vertex_id> > cid_edges;
    std::vector< std::pair<vertex_id, vertex_id> > cid2_edges;

    for(unsigned edge=0; edge<3; ++edge)
    {
      cid_edges.push_back(std::make_pair(cid_vertices[edge], cid_vertices[(edge+1)%3]));
      cid2_edges.push_back(std::make_pair(cid2_vertices[edge], cid2_vertices[(edge+1)%3]));
    }

    // Map vertices and midpoints of cid2 onto degrees of freedom
    unsigned cid2_dof=0;

    std::map<vertex_id, unsigned> cid2_vertex_dof;
    for(unsigned v=0; v<cid2_vertices.size(); ++v)
    {
      cid2_vertex_dof[cid2_vertices[v]] = cid2_dof;
      ++cid2_dof;
    }

    std::map<std::pair<vertex_id, vertex_id>, unsigned, unordered_pair_compare<vertex_id> > cid2_midpoint_dof;
    for(unsigned m=0; m<cid_edges.size(); ++m)
    {
      cid2_midpoint_dof[cid2_edges[m]] = cid2_dof;
      ++cid2_dof;
    }

    // Iterate over degrees of freedom on cid and find ones that correspond to common vertices and midpoints
    std::vector< std::pair<unsigned, unsigned> > common;
    unsigned cid_dof = 0;

    for(unsigned v=0; v<cid_vertices.size(); ++v)
    {
      const std::map<vertex_id, unsigned>::const_iterator sharedVertexIter = cid2_vertex_dof.find(cid_vertices[v]);

      if (sharedVertexIter != cid2_vertex_dof.end())
        common.push_back(std::make_pair(cid_dof, sharedVertexIter->second));

      ++cid_dof;
    }

    for(unsigned m=0; m<cid_edges.size(); ++m)
    {
      const std::map<std::pair<vertex_id, vertex_id>, unsigned>::const_iterator sharedEdgeIter = cid2_midpoint_dof.find(cid_edges[m]);

      if (sharedEdgeIter != cid2_midpoint_dof.end())
      {
        common.push_back(std::make_pair(cid_dof, sharedEdgeIter->second));
      }

      ++cid_dof;
    }
    return common;
  }

};

}

#endif
