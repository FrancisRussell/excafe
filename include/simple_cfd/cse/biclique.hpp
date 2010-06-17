#ifndef SIMPLE_CFD_CSE_BICLIQUE_HPP
#define SIMPLE_CFD_CSE_BICLIQUE_HPP

#include <set>
#include <utility>
#include <boost/foreach.hpp>
#include "properties.hpp"
#include <simple_cfd/util/maybe.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace cse
{

template<typename G>
class Biclique
{
public:
  typedef G graph_t;
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<graph_t>::edge_descriptor   edge_descriptor;

protected:
  friend class BicliqueSearchSpace< Biclique<graph_t> >;

  graph_t* graph;
  std::set<vertex_descriptor> cubeVertices;
  std::set<vertex_descriptor> coKernelVertices;
  int cubeValueSum;
  int coKernelValueSum;

  class EdgeRemover
  {
  private:
    const graph_t* graph;
    std::set< std::pair<std::size_t, std::size_t> > termIDs;

  public:
    EdgeRemover(const graph_t& _graph, const std::set< std::pair<std::size_t, std::size_t> > _termIDs) : 
      graph(&_graph), termIDs(_termIDs)
    {
    }

    bool operator()(const edge_descriptor& e) const
    {
      return termIDs.find(get(term_id(), *graph, e)) != termIDs.end();
    }
  };

  template<typename InputIterator>
  static int getValue(const graph_t& graph, const InputIterator begin, const InputIterator end)
  {
    int result = 0;

    BOOST_FOREACH(const vertex_descriptor& v, std::make_pair(begin, end))
    {
      result += get(mul_count(), graph, v);
    }

    return result;
  }

  void removeUnconnected(const vertex_descriptor& v, std::set<vertex_descriptor>& vertices)
  {
    std::set<vertex_descriptor> newVertices;

    BOOST_FOREACH(const edge_descriptor& edge, out_edges(v, *graph))
    {
      const typename std::set<vertex_descriptor>::const_iterator vIter = vertices.find(target(edge, *graph));
      if (vIter != vertices.end())
      {
        newVertices.insert(*vIter);
      }
    }

    std::swap(vertices, newVertices);
  }

  static int getValue(const int cubeValueSum, const int numCubes, const int coKernelValueSum, const int numCoKernels)
  {
    const int multiplyWeight = 1;

    if (numCubes == 0 || numCoKernels == 0)
      return 0;

    const int origAdds = (numCubes - 1) * numCoKernels;
    const int origMuls = cubeValueSum * numCoKernels + coKernelValueSum * numCubes;

    const int rectAdds = numCubes - 1;
    const int rectMuls = cubeValueSum + coKernelValueSum;

    const int savedAdds = origAdds - rectAdds;
    const int savedMuls = origMuls - rectMuls;

    return multiplyWeight*savedMuls + savedAdds;
  }

  Biclique(const Biclique& parent, const vertex_descriptor v) : graph(parent.graph),
    cubeVertices(parent.cubeVertices), coKernelVertices(parent.coKernelVertices)
  {
    const bool isCube = get(is_cube(), *graph, v);
    const bool firstVertex = isCube ? cubeVertices.empty() : coKernelVertices.empty();

    if (isCube)
    {
      cubeVertices.insert(v);
      removeUnconnected(v, coKernelVertices);
    }
    else
    {
      coKernelVertices.insert(v);
      removeUnconnected(v, cubeVertices);
    }

    if (firstVertex)
    {
      BOOST_FOREACH(const edge_descriptor& edge, out_edges(v, *graph))
      {
        (isCube ? coKernelVertices : cubeVertices).insert(target(edge, *graph));
      }
    }

    calculateValue();
  }

  void calculateValue()
  {
    cubeValueSum = getValue(*graph, cubeVertices.begin(), cubeVertices.end());
    coKernelValueSum = getValue(*graph, coKernelVertices.begin(), coKernelVertices.end());
  }

public:
  Biclique(graph_t& _graph) : graph(&_graph), cubeValueSum(0), coKernelValueSum(0)
  {
  }

  const graph_t& getGraph() const
  {
    return *graph;
  }

  void print() const
  {
    std::cout << "num_cubes=" << cubeVertices.size() << ", num_cokernels=" << coKernelVertices.size();
    std::cout << ", value=" << getValue() << std::endl;
  }

  bool empty() const
  {
    return cubeVertices.empty() && coKernelVertices.empty();
  }

  int getValue() const
  {
    return getValue(cubeValueSum, cubeVertices.size(), coKernelValueSum, coKernelVertices.size());
  }

  std::size_t numCubes() const
  {
    return cubeVertices.size();
  }

  std::size_t numCoKernels() const
  {
    return coKernelVertices.size();
  }

  int getCubeValueSum() const
  {
    return cubeValueSum;
  }

  int getCoKernelValueSum() const
  {
    return coKernelValueSum;
  }

  std::set<vertex_descriptor> getCubeVertices() const
  {
    return cubeVertices;
  }

  std::set<vertex_descriptor> getCoKernelVertices() const
  {
    return coKernelVertices;
  }

  Biclique addVertex(const vertex_descriptor& v) const
  {
    return Biclique(*this, v);
  }

  void swap(Biclique& b)
  {
    std::swap(graph, b.graph);
    std::swap(cubeVertices, b.cubeVertices);
    std::swap(coKernelVertices, b.coKernelVertices);
  }

  SOP getSOP() const
  {
    SOP result;
    BOOST_FOREACH(const vertex_descriptor& v, cubeVertices)
    {
      result.append(get(term_cube(), *graph, v));
    }
    return result;
  }

  void rewriteAndUpdate(const vertex_descriptor& newCubeVertex, SOP* const sops)
  {
    const Cube newCube = get(term_cube(), *graph, newCubeVertex);
    std::set<std::pair<std::size_t, std::size_t> > termIDs;

    std::cout << "performing term update" << std::endl;
    BOOST_FOREACH(const vertex_descriptor& coKernelVertex, coKernelVertices)
    {
      // TODO: get polynomialID from coKernelVertex instead of assigning it multiple times in the loop.
      std::size_t polynomialID = 0;
      BOOST_FOREACH(const edge_descriptor& edge, out_edges(coKernelVertex, *graph))
      {
        const vertex_descriptor cubeVertex = target(edge, *graph);
        if (cubeVertices.find(cubeVertex) != cubeVertices.end())
        {
          const std::pair<std::size_t, std::size_t> termID = get(term_id(), *graph, edge);
          polynomialID = termID.first;
          
          const bool inserted = termIDs.insert(termID).second;
          if (!inserted)
            CFD_EXCEPTION("Duplicate term found in biclique.");

          const bool deleted = sops[polynomialID].deleteTerm(termID.second);
          if (!deleted)
            CFD_EXCEPTION("Failed to remove factorised term from SOP.");
        }
      }

      // Add new term to polynomial
      const std::size_t termID = sops[polynomialID].append(get(term_cokernel(), *graph, coKernelVertex) + newCube);

      // Add new edge to graph and tag with polynomial and termID;
      const std::pair<edge_descriptor, bool> edgePair = add_edge(coKernelVertex, newCubeVertex, *graph);
      if (!edgePair.second)
        CFD_EXCEPTION("Tried to insert duplicate edge when updating KCM. This should never happen.");
      put(term_id(), *graph, edgePair.first, std::make_pair(polynomialID, termID));

    }

    // Remove old terms from graph
    remove_edge_if(EdgeRemover(*graph, termIDs), *graph);
  }
};

}

}

namespace std
{

template<typename G>
void swap(cfd::cse::Biclique<G>& a, cfd::cse::Biclique<G>& b)
{
  a.swap(b);
}

}

#endif
