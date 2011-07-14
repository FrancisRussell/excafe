#ifndef SIMPLE_CFD_CSE_KCM_HPP
#define SIMPLE_CFD_CSE_KCM_HPP

// Avoid Boost Graph including deprecated hash_set header
#define BOOST_NO_HASH

#include <cstddef>
#include <vector>
#include <map>
#include <utility>
#include <queue>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "cse_fwd.hpp"
#include "cube.hpp"
#include "sop.hpp"
#include "properties.hpp"
#include "polynomial_index.hpp"
#include "sop_map.hpp"

namespace cfd
{

namespace cse
{

class KCM
{
private:
  typedef boost::property<term_cube, Cube,
          boost::property<term_cokernel, Cube,
          boost::property<mul_count, int,
          boost::property<is_cube, bool,
          boost::property<is_unit, bool,
          boost::property<is_numeric, bool,
          boost::property<polynomial_id, PolynomialIndex,
          boost::property<cube_ordering, std::pair<int, unsigned>
          > > > > > > > > VertexProperty;

  // std::pair<polynomial_id, term_number>
  typedef boost::property< term_id, std::pair<PolynomialIndex, std::size_t> > EdgeProperty;

  typedef boost::adjacency_list<
    boost::vecS,
    boost::listS,
    boost::undirectedS,
    VertexProperty,
    EdgeProperty,
    boost::no_property,
    boost::listS> graph_t;

  typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<graph_t>::edge_descriptor   edge_descriptor;
  typedef boost::graph_traits<graph_t>::out_edge_iterator out_edge_iterator;
  typedef Biclique<graph_t> biclique_t;
  typedef BicliqueSearch<graph_t> biclique_search_t;

  NewLiteralCreator& literalCreator;
  SOPMap& sops;
  std::map<Cube, vertex_descriptor> cubeVertices;
  graph_t graph;

  template<typename PriorityQueue>
  void addSearchSpaces(PriorityQueue& out)
  {
    std::size_t count=0;
    BOOST_FOREACH(const vertex_descriptor v, vertices(graph))
    {
      if (get(is_cube(), graph, v))
      {
        // out_degree is potentially O(n) so we compare iterators instead.
        out_edge_iterator edgesBegin, edgesEnd;
        boost::tie(edgesBegin, edgesEnd) = out_edges(v, graph);

        if (edgesBegin != edgesEnd)
        {
          out.push(biclique_search_t(graph, v));
          ++count;
        }
      }
    }
    std::cout << "Added " << count << " cubes to search space." << std::endl;
  }

  void orderCubes();
  PolynomialIndex addPolynomial(const SOP& sop);
  void addPolynomial(const PolynomialIndex& polynomialID);

public:
  typedef std::vector<SOP>::const_iterator iterator;
  typedef std::vector<SOP>::const_iterator const_iterator;

  KCM(NewLiteralCreator& _literalCreator);

  vertex_descriptor addCoKernel(const PolynomialIndex& polynomialID, const Cube& cokernel);
  vertex_descriptor addCube(const Cube& c);
  bool factorise();
  std::size_t numEdges() const;
  std::size_t numCubes() const;
  std::size_t numCoKernels() const;
  std::size_t numAdditions() const;
  std::size_t numMultiplies() const;
  void updateGraph(const Biclique<graph_t>& biclique);
  void removeBiclique(const Biclique<graph_t>& biclique);
};

}

}

#endif
