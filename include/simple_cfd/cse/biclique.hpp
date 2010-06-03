#ifndef SIMPLE_CFD_CSE_BICLIQUE_HPP
#define SIMPLE_CFD_CSE_BICLIQUE_HPP

#include <set>
#include <boost/foreach.hpp>
#include "properties.hpp"

namespace cfd
{

namespace cse
{

template<typename G>
class Biclique
{
private:
  typedef G graph_t;
  typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<graph_t>::edge_descriptor   edge_descriptor;

  graph_t* graph;
  std::set<vertex_descriptor> cubeVertices;
  std::set<vertex_descriptor> coKernelVertices;

  int getValue(const std::set<vertex_descriptor>& vertices) const
  {
    int result = 0;

    BOOST_FOREACH(const vertex_descriptor& v, vertices)
    {
      result += get(mul_count(), *graph, v);
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

public:
  Biclique(graph_t& _graph) : graph(&_graph)
  {
  }

  Biclique(graph_t& _graph, const vertex_descriptor seed) : graph(&_graph)
  {
    addVertex(seed);

    BOOST_FOREACH(const edge_descriptor& edge, out_edges(seed, *graph))
    {
      addVertex(target(edge, *graph));
    }
  }

  int getValue() const
  {
    const int multiplyWeight = 1;
    const int numCubes = cubeVertices.size();
    const int numCoKernels = coKernelVertices.size();

    if (numCubes == 0 || numCoKernels == 0)
      return 0;

    const int cubeValueSum = getValue(cubeVertices);
    const int coKernelValueSum = getValue(coKernelVertices);

    const int origAdds = (numCubes - 1) * numCoKernels;
    const int origMuls = cubeValueSum * numCoKernels + coKernelValueSum * numCubes;

    const int rectAdds = numCubes - 1;
    const int rectMuls = cubeValueSum + coKernelValueSum;

    const int savedAdds = origAdds - rectAdds;
    const int savedMuls = origMuls - rectMuls;

    return multiplyWeight*savedMuls + savedAdds;
  }

  void addVertex(const vertex_descriptor& v)
  {
    const bool isCube = get(is_cube(), *graph, v);

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
