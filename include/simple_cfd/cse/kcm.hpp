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
#include <boost/graph/adjacency_list.hpp>
#include "cube.hpp"
#include "sop.hpp"
#include "properties.hpp"
#include "biclique.hpp"
#include "biclique_search_space.hpp"
#include "polynomial_index.hpp"
#include <simple_cfd/exception.hpp>

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
          boost::property<is_cube, bool
          > > > > VertexProperty;

  // std::pair<polynomial_id, term_number>
  typedef boost::property< term_id, std::pair<std::size_t, std::size_t> > EdgeProperty;

  typedef boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    VertexProperty,
    EdgeProperty,
    boost::no_property,
    boost::listS> graph_t;

  typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<graph_t>::edge_descriptor   edge_descriptor;
  typedef Biclique<graph_t> biclique_t;
  typedef BicliqueSearchSpace<biclique_t> biclique_search_t;

  NewLiteralCreator& literalCreator;
  std::vector<SOP> polynomials;
  std::map<Cube, vertex_descriptor> cubeVertices;
  graph_t graph;

  BicliqueSearchSpace<biclique_t> getSearchSpace(const vertex_descriptor v)
  {
    std::set<vertex_descriptor> adjacent;
    BOOST_FOREACH(const edge_descriptor& edge, out_edges(v, graph))
    {
      adjacent.insert(target(edge, graph));
    }

    std::set<vertex_descriptor> neighbours;
    BOOST_FOREACH(const vertex_descriptor& a, adjacent)
    {
      BOOST_FOREACH(const edge_descriptor& edge, out_edges(a, graph))
      {
        neighbours.insert(target(edge, graph));
      }
    }

    const biclique_t seed(biclique_t(graph).addVertex(v));
    biclique_search_t result(seed, neighbours.begin(), neighbours.end());
    return result;
  }

  template<typename PriorityQueue>
  void addSearchSpaces(PriorityQueue& out)
  {
    BOOST_FOREACH(const vertex_descriptor v, vertices(graph))
    {
      if (get(is_cube(), graph, v))
      {
        std::cout << "Adding search space for vertex " << v << "..." << std::endl;
        out.push(getSearchSpace(v));
      }
    }
  }

public:
  typedef std::vector<SOP>::const_iterator iterator;
  typedef std::vector<SOP>::const_iterator const_iterator;

  const_iterator begin() const
  {
    return polynomials.begin();
  }

  const_iterator end() const
  {
    return polynomials.end();
  }

  KCM(NewLiteralCreator& _literalCreator) : literalCreator(_literalCreator)
  {
  }

  std::size_t addPolynomial(const SOP& sop)
  {
    const std::size_t polynomialID = polynomials.size();
    polynomials.push_back(sop);

    const SOP::kernel_set_t kernels = sop.getKernels();
    BOOST_FOREACH(const SOP::kernel_set_t::value_type kernel, kernels)
    {
      const vertex_descriptor coKernelVertex = addCoKernel(kernel.second);
      const SOP& sop = kernel.first;

      for(SOP::const_iterator iter = sop.begin(); iter != sop.end(); ++iter)
      {
        const vertex_descriptor cubeVertex = addCube(*iter);
        const std::size_t termID = sop.getTermNumber(iter);
        const std::pair<edge_descriptor, bool> edgePair = add_edge(coKernelVertex, cubeVertex, graph);

        if (!edgePair.second)
          CFD_EXCEPTION("Attemped to insert duplicate edge into KCM. This should never happen.");

        const edge_descriptor edge = edgePair.first;
        put(term_id(), graph, edge, std::make_pair(polynomialID, termID));
      }

      std::cout << "kernel: " << kernel.first << ", co-kernel: " << kernel.second << std::endl;
    }

    return polynomialID;
  }

  vertex_descriptor addCoKernel(const Cube& cokernel)
  {
    const vertex_descriptor v = add_vertex(graph);
    put(is_cube(), graph, v, false);
    put(term_cokernel(), graph, v, cokernel);
    put(mul_count(), graph, v, cokernel.numMultiplies());
    return v;
  }

  vertex_descriptor addCube(const Cube& c)
  {
    const std::map<Cube, vertex_descriptor>::const_iterator vertexIter = cubeVertices.find(c);

    if (vertexIter != cubeVertices.end())
    {
      return vertexIter->second;
    }
    else
    {
      const vertex_descriptor v = add_vertex(graph);
      cubeVertices.insert(std::make_pair(c, v));
      put(is_cube(), graph, v, true);
      put(term_cube(), graph, v, c);
      put(mul_count(), graph, v, c.numMultiplies());
      return v;
    }
  }

  void factorise()
  {
    typedef BicliqueSearchSpaceComparator<biclique_t> biclique_search_comparator_t;
    std::priority_queue<biclique_search_t, std::vector<biclique_search_t>, biclique_search_comparator_t> queue;
    addSearchSpaces(queue);
    biclique_t best(graph);

    while(!queue.empty())
    {
      const biclique_search_t bs = queue.top();
      queue.pop();

      bs.print();

      if (!bs.finished() && bs.getMaximalValue() > best.getValue())
      {
        const std::pair<biclique_search_t, biclique_search_t> pair = bs.split();
        queue.push(pair.first);
        queue.push(pair.second);
      }

      if (best.getValue() < bs.getValue())
      {
        best = bs.getBiclique();
        std::cout << "New best score: " << best.getValue() << std::endl;
        bs.print();
      }
    }

    removeBiclique(best);
  }

  std::size_t numEdges() const
  {
    return num_edges(graph);
  }

  std::size_t numCubes() const
  {
    std::size_t count = 0;
    BOOST_FOREACH(const vertex_descriptor& v, vertices(graph))
    {
      if (get(is_cube(), graph, v))
      {
        ++count;
      }
    }
    return count;
  }

  std::size_t numCoKernels() const
  {
    return num_vertices(graph) - numCubes();
  }

  void removeBiclique(const Biclique<graph_t>& biclique)
  {
    //FIXME: check for duplicate terms in biclique.

    const SOP newSOP = biclique.getSOP();
    const std::size_t newSOPIndex = addPolynomial(newSOP);
    const unsigned literal = literalCreator.getLiteralID(PolynomialIndex(newSOPIndex));
    const std::map<std::size_t, SOPRewrite> rewrites = biclique.getRewrites(literal);

    typedef std::pair<std::size_t, SOPRewrite> rewrite_mapping_t;
    BOOST_FOREACH(const rewrite_mapping_t& rewrite, rewrites)
    {
      polynomials[rewrite.first] = rewrite.second(polynomials[rewrite.first]);
    }
  }
};

}

}

#endif
