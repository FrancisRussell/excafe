#ifndef EXCAFE_CSE_BICLIQUE_HPP
#define EXCAFE_CSE_BICLIQUE_HPP

#include <set>
#include <utility>
#include <boost/foreach.hpp>
#include <boost/operators.hpp>
#include <boost/container/flat_set.hpp>
#include "properties.hpp"
#include "sop_map.hpp"
#include "vertex_info.hpp"
#include "polynomial_index.hpp"
#include <excafe/util/lazy_copy.hpp>
#include <excafe/exception.hpp>

namespace excafe
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
  graph_t* graph;
  util::LazyCopy< boost::container::flat_set<vertex_descriptor> > cubeVertices;
  util::LazyCopy< boost::container::flat_set<vertex_descriptor> > coKernelVertices;
  VertexInfo<graph_t> cubeInfo;
  VertexInfo<graph_t> coKernelInfo;

  template<typename InputIterator>
  static VertexInfo<graph_t> getValue(const graph_t& graph, const InputIterator begin, const InputIterator end)
  {
    VertexInfo<graph_t> result;
    BOOST_FOREACH(const vertex_descriptor& v, std::make_pair(begin, end))
      result.addVertex(graph, v);

    return result;
  }

  void removeUnconnected(const vertex_descriptor& v, boost::container::flat_set<vertex_descriptor>& vertices)
  {
    boost::container::flat_set<vertex_descriptor> newVertices;

    BOOST_FOREACH(const edge_descriptor& edge, out_edges(v, *graph))
    {
      const typename boost::container::flat_set<vertex_descriptor>::const_iterator vIter = vertices.find(target(edge, *graph));
      if (vIter != vertices.end())
        newVertices.insert(*vIter);
    }

    std::swap(vertices, newVertices);
  }

  static int getValue(const VertexInfo<graph_t>& cubeInfo, const VertexInfo<graph_t>& coKernelInfo)
  {
    const int multiplyWeight = 1;

    if (cubeInfo.num() == 0 || coKernelInfo.num() == 0)
      return 0;

    const int origAdds = (cubeInfo.num() - 1) * coKernelInfo.num();
    const int origMuls = cubeInfo.getValue() * coKernelInfo.num() +
                         coKernelInfo.getValue() * cubeInfo.num() +
                         coKernelInfo.numNonUnit() * cubeInfo.numNonUnit() -
                         coKernelInfo.numHaveCoefficients() * cubeInfo.numHaveCoefficients();

    const int rectAdds = cubeInfo.num() - 1;
    const int rectMuls = cubeInfo.getValue() + 
                         coKernelInfo.getValue() + 
                         coKernelInfo.numNonUnit();

    const int savedAdds = origAdds - rectAdds;
    const int savedMuls = origMuls - rectMuls;

    return multiplyWeight*savedMuls + savedAdds;
  }

  Biclique(const Biclique& parent, const vertex_descriptor v) : graph(parent.graph),
    cubeVertices(parent.cubeVertices), coKernelVertices(parent.coKernelVertices)
  {
    const bool isCube = get(is_cube(), *graph, v);
    const bool firstVertex = isCube ? cubeVertices->empty() : coKernelVertices->empty();

    if (isCube)
    {
      cubeVertices->insert(v);
      removeUnconnected(v, *coKernelVertices);
    }
    else
    {
      coKernelVertices->insert(v);
      removeUnconnected(v, *cubeVertices);
    }

    if (firstVertex)
    {
      BOOST_FOREACH(const edge_descriptor& edge, out_edges(v, *graph))
        (isCube ? coKernelVertices : cubeVertices)->insert(target(edge, *graph));
    }

    calculateValue();
  }

  void calculateValue()
  {
    cubeInfo = getValue(*graph, cubeVertices->begin(), cubeVertices->end());
    coKernelInfo = getValue(*graph, coKernelVertices->begin(), coKernelVertices->end());
  }

public:
  Biclique(graph_t& _graph) : graph(&_graph)
  {
  }

  const graph_t& getGraph() const
  {
    return *graph;
  }

  void print() const
  {
    std::cout << "num_cubes=" << numCubes() << ", num_cokernels=" << numCoKernels();
    std::cout << ", value=" << getValue() << std::endl;
  }

  bool empty() const
  {
    return cubeVertices->empty() && coKernelVertices->empty();
  }

  int getValue() const
  {
    return getValue(cubeInfo, coKernelInfo);
  }

  std::size_t numCubes() const
  {
    return cubeVertices->size();
  }

  std::size_t numCoKernels() const
  {
    return coKernelVertices->size();
  }

  int getCubeValueSum() const
  {
    return cubeInfo.getValue();
  }

  int getCoKernelValueSum() const
  {
    return coKernelInfo.getValue();
  }

  std::set<vertex_descriptor> getCubeVertices() const
  {
    return std::set<vertex_descriptor>(cubeVertices->begin(), cubeVertices->end());
  }

  std::set<vertex_descriptor> getCoKernelVertices() const
  {
    return std::set<vertex_descriptor>(coKernelVertices->begin(), coKernelVertices->end());
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
    std::swap(cubeInfo, b.cubeInfo);
    std::swap(coKernelInfo, b.coKernelInfo);
  }

  SOP getSOP() const
  {
    SOP result;
    BOOST_FOREACH(const vertex_descriptor& v, *cubeVertices)
      result.append(get(term_cube(), *graph, v));
    return result;
  }

  std::set<PolynomialIndex> getModifiedPolynomials() const
  {
    std::set<PolynomialIndex> polynomials;
    BOOST_FOREACH(const vertex_descriptor& coKernelVertex, *coKernelVertices)
    {
      const PolynomialIndex polynomialID = get(polynomial_id(), *graph, coKernelVertex);
      polynomials.insert(polynomialID);
    }
    return polynomials;
  }

  void rewritePolynomials(const vertex_descriptor& newCubeVertex, SOPMap& sops) const
  {
    const Cube newCube = get(term_cube(), *graph, newCubeVertex);
    std::set<std::pair<PolynomialIndex, std::size_t> > termIDs;

    BOOST_FOREACH(const vertex_descriptor& coKernelVertex, *coKernelVertices)
    {
      const PolynomialIndex polynomialID = get(polynomial_id(), *graph, coKernelVertex);
      BOOST_FOREACH(const edge_descriptor& edge, out_edges(coKernelVertex, *graph))
      {
        const vertex_descriptor cubeVertex = target(edge, *graph);
        if (cubeVertices->find(cubeVertex) != cubeVertices->end())
        {
          const std::size_t termID = get(term_id(), *graph, edge);

          const bool inserted = termIDs.insert(std::make_pair(polynomialID, termID)).second;
          if (!inserted)
            CFD_EXCEPTION("Duplicate term found in biclique.");

          const bool deleted = sops[polynomialID].deleteTerm(termID);
          if (!deleted)
            CFD_EXCEPTION("Failed to remove factorised term from SOP.");
        }
      }

      // Add new term to polynomial
      sops[polynomialID].append(get(term_cube(), *graph, coKernelVertex) + newCube);
    }
  }
};

}

}

namespace std
{

template<typename G>
void swap(excafe::cse::Biclique<G>& a, excafe::cse::Biclique<G>& b)
{
  a.swap(b);
}

}

#endif
