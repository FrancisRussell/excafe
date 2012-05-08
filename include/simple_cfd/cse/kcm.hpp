#ifndef SIMPLE_CFD_CSE_KCM_HPP
#define SIMPLE_CFD_CSE_KCM_HPP

// Avoid Boost Graph including deprecated hash_set header
#define BOOST_NO_HASH

#include <cstddef>
#include <cassert>
#include <vector>
#include <map>
#include <set>
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
#include "biclique_search.hpp"

namespace cfd
{

namespace cse
{

class KCM
{
private:
  typedef boost::property<term_cube, Cube,
          boost::property<mul_count, int,
          boost::property<is_cube, bool,
          boost::property<is_unit, bool,
          boost::property<is_numeric, bool,
          boost::property<has_coefficient, bool,
          boost::property<polynomial_id, PolynomialIndex,
          boost::property<cube_ordering, std::pair<int, unsigned>
          > > > > > > > > VertexProperty;

  // term_number
  typedef boost::property<term_id, std::size_t> EdgeProperty;

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

  class CubeComparator
  {
  private:
    const graph_t* graph;

  public:
    CubeComparator(const graph_t& _graph) : graph(&_graph)
    {
    }

    bool operator()(const vertex_descriptor& a, const vertex_descriptor& b) const
    {
      return get(cube_ordering(), *graph, a) < get(cube_ordering(), *graph, b);
    }
  };

  NewLiteralCreator& literalCreator;
  SOPMap& sops;
  std::map<Cube, vertex_descriptor> cubeVertices;
  graph_t graph;

  template<typename PriorityQueue>
  void addSearchSpaces(PriorityQueue& out)
  {
    // Construct set of cubes ordered by the same property we use to
    // choose which cube to next grow a biclique by.
    const CubeComparator cubeComparator(graph);
    std::set<vertex_descriptor, CubeComparator> orderedCubes(cubeComparator);

    BOOST_FOREACH(const vertex_descriptor v, vertices(graph))
    {
      if (get(is_cube(), graph, v))
      {
        const bool inserted = orderedCubes.insert(v).second;
        assert(inserted && "Cube ordering value needs to be unique");
      }
    }

    // Add search spaces to priority queue only if their set of
    // co-kernels have not been seen before. If they have, the previous
    // search space contains the one just constructed.
    std::set< std::set<vertex_descriptor> > coKernelSets;
    std::size_t count=0;
    BOOST_FOREACH(const vertex_descriptor v, orderedCubes)
    {
      // out_degree is potentially O(n) so we compare iterators instead.
      out_edge_iterator edgesBegin, edgesEnd;
      boost::tie(edgesBegin, edgesEnd) = out_edges(v, graph);

      if (edgesBegin != edgesEnd)
      {
        const biclique_search_t searchSpace(graph, v);
        const bool isNewSearchSpace = coKernelSets.insert(searchSpace.getCoKernelVertices()).second;
        if (isNewSearchSpace)
        {
          out.push(searchSpace);
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
