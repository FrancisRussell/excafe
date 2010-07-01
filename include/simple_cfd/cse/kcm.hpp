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
#include <boost/utility.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "cube.hpp"
#include "sop.hpp"
#include "properties.hpp"
#include "biclique.hpp"
#include "biclique_search.hpp"
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
          boost::property<is_cube, bool,
          boost::property<is_one, bool,
          boost::property<polynomial_id, std::size_t,
          boost::property<cube_ordering, std::pair<int, unsigned>
          > > > > > > > VertexProperty;

  // std::pair<polynomial_id, term_number>
  typedef boost::property< term_id, std::pair<std::size_t, std::size_t> > EdgeProperty;

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
  std::vector<SOP> polynomials;
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

  void orderCubes()
  {
    unsigned id=0;
    BOOST_FOREACH(const vertex_descriptor& vertex, vertices(graph))
    {
      if (get(is_cube(), graph, vertex))
      {
        put(cube_ordering(), graph, vertex, std::make_pair(out_degree(vertex, graph), id));
        ++id;
      }
    }
  }

  void addPolynomial(const std::size_t polynomialID)
  {
    const SOP& sop = polynomials[polynomialID];
    const SOP::kernel_set_t kernels = sop.getKernels();
    BOOST_FOREACH(const SOP::kernel_set_t::value_type kernel, kernels)
    {
      const vertex_descriptor coKernelVertex = addCoKernel(polynomialID, kernel.second);
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
      //std::cout << "kernel: " << kernel.first << ", co-kernel: " << kernel.second << std::endl;
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
    addPolynomial(polynomialID);
    return polynomialID;
  }

  vertex_descriptor addCoKernel(const std::size_t polynomialID, const Cube& cokernel)
  {
    const vertex_descriptor v = add_vertex(graph);
    put(is_cube(), graph, v, false);
    put(polynomial_id(), graph, v, polynomialID);
    put(term_cokernel(), graph, v, cokernel);
    put(mul_count(), graph, v, cokernel.numMultiplies());
    put(is_one(), graph, v, cokernel.isOne());
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
      put(is_one(), graph, v, c.isOne());
      return v;
    }
  }

  bool factorise()
  {
    orderCubes();

    std::priority_queue<biclique_search_t, std::vector<biclique_search_t>, BicliqueSearchComparator> queue;
    addSearchSpaces(queue);
    biclique_t best(graph);

    while(!queue.empty())
    {
      const biclique_search_t bs = queue.top();
      queue.pop();

      if (bs.getMaximalValue() <= best.getValue())
      {
        break;
      }
      else
      {
        //bs.print();
        //std::cout << std::endl;

        if (!bs.isFinished())
        {
          const std::pair<biclique_search_t, biclique_search_t> pair = bs.split();

          if (pair.first.getMaximalValue() > best.getValue())
            queue.push(pair.first);

          if (pair.second.getMaximalValue() > best.getValue())
            queue.push(pair.second);
        }

        if (best.getValue() < bs.getValue())
        {
          best = bs.getBiclique();
          std::cout << "New best score: " << best.getValue() << std::endl;
        }
      }
    }

    if (best.getValue() > 0)
    {
      removeBiclique(best);
      return true;
    }
    else
    {
      return false;
    }
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
        ++count;
    }
    return count;
  }

  std::size_t numCoKernels() const
  {
    return num_vertices(graph) - numCubes();
  }

  std::size_t numAdditions() const
  {
    std::size_t result=0;
    BOOST_FOREACH(const SOP& sop, polynomials)
    {
      result += sop.numAdditions();
    }
    return result;
  }

  std::size_t numMultiplies() const
  {
    std::size_t result=0;
    BOOST_FOREACH(const SOP& sop, polynomials)
    {
      result += sop.numMultiplies();
    }
    return result;
  }

  void updateGraph(const Biclique<graph_t>& biclique)
  {
    const std::set<std::size_t> modifiedPolynomials(biclique.getModifiedPolynomials());

    typedef boost::graph_traits<graph_t>::vertex_iterator vertex_iter;
    vertex_iter vi, viEnd;

    boost::tie(vi, viEnd) = vertices(graph);
    while(vi != viEnd)
    {
      const vertex_iter viNext = boost::next(vi);
      if (!get(is_cube(), graph, *vi))
      {
        if (modifiedPolynomials.find(get(polynomial_id(), graph, *vi)) != modifiedPolynomials.end())
        {
          clear_vertex(*vi, graph);
          remove_vertex(*vi, graph);
        }
      }

      vi = viNext;
    }

    BOOST_FOREACH(const std::size_t polynomialID, modifiedPolynomials)
    {
      addPolynomial(polynomialID);
    }
  }

  void removeBiclique(const Biclique<graph_t>& biclique)
  {
    const SOP newSOP = biclique.getSOP();
    const std::size_t newSOPIndex = addPolynomial(newSOP);
    const unsigned literal = literalCreator.getLiteralID(PolynomialIndex(newSOPIndex));

    // Rewrite polynomials
    biclique.rewritePolynomials(addCube(literal), &polynomials[0]);

    // Update graph
    updateGraph(biclique);
  }
};

}

}

#endif
