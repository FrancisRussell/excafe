#ifndef SIMPLE_CFD_CSE_BICLIQUE_HPP
#define SIMPLE_CFD_CSE_BICLIQUE_HPP

#include <set>
#include <utility>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include "properties.hpp"
#include "sop_map.hpp"
#include "polynomial_index.hpp"
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
  graph_t* graph;
  std::set<vertex_descriptor> cubeVertices;
  std::set<vertex_descriptor> coKernelVertices;
  std::size_t nonOneCubes;
  std::size_t nonOneCoKernels;
  int cubeValueSum;
  int coKernelValueSum;

  template<typename InputIterator>
  static std::pair<std::size_t, int> getValue(const graph_t& graph, const InputIterator begin, const InputIterator end)
  {
    std::size_t nonOneCount = 0;
    int value = 0;

    BOOST_FOREACH(const vertex_descriptor& v, std::make_pair(begin, end))
    {
      if (!get(is_unit(), graph, v)) 
        ++nonOneCount;

      value += get(mul_count(), graph, v);
    }

    return std::make_pair(nonOneCount, value);
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

  static int getValue(const int cubeValueSum, const std::size_t numNonOneCubes, const std::size_t numCubes, 
    const int coKernelValueSum, const std::size_t numNonOneCoKernels, const std::size_t numCoKernels)
  {
    const int multiplyWeight = 1;

    if (numCubes == 0 || numCoKernels == 0)
      return 0;

    const int origAdds = (numCubes - 1) * numCoKernels;
    const int origMuls = cubeValueSum * numCoKernels + coKernelValueSum * numCubes + numNonOneCubes*numNonOneCoKernels;

    const int rectAdds = numCubes - 1;
    const int rectMuls = cubeValueSum + coKernelValueSum + numNonOneCoKernels;

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
    boost::tie(nonOneCubes, cubeValueSum) = getValue(*graph, cubeVertices.begin(), cubeVertices.end());
    boost::tie(nonOneCoKernels, coKernelValueSum) = getValue(*graph, coKernelVertices.begin(), coKernelVertices.end());
  }

public:
  Biclique(graph_t& _graph) : graph(&_graph), nonOneCubes(0), nonOneCoKernels(0), 
    cubeValueSum(0), coKernelValueSum(0)
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
    return getValue(cubeValueSum, nonOneCubes, cubeVertices.size(), 
                    coKernelValueSum, nonOneCoKernels, coKernelVertices.size());
  }

  std::size_t numCubes() const
  {
    return cubeVertices.size();
  }

  std::size_t numNonOneCubes() const
  {
    return nonOneCubes;
  }

  std::size_t numCoKernels() const
  {
    return coKernelVertices.size();
  }

  std::size_t numNonOneCoKernels() const
  {
    return nonOneCoKernels;
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
    std::swap(cubeValueSum, b.cubeValueSum);
    std::swap(coKernelValueSum, b.coKernelValueSum);
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

  std::set<PolynomialIndex> getModifiedPolynomials() const
  {
    std::set<PolynomialIndex> polynomials;
    BOOST_FOREACH(const vertex_descriptor& coKernelVertex, coKernelVertices)
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

    BOOST_FOREACH(const vertex_descriptor& coKernelVertex, coKernelVertices)
    {
      const PolynomialIndex polynomialID = get(polynomial_id(), *graph, coKernelVertex);
      BOOST_FOREACH(const edge_descriptor& edge, out_edges(coKernelVertex, *graph))
      {
        const vertex_descriptor cubeVertex = target(edge, *graph);
        if (cubeVertices.find(cubeVertex) != cubeVertices.end())
        {
          const std::pair<PolynomialIndex, std::size_t> termID = get(term_id(), *graph, edge);

          // If the polynomial ID of the term associated with this edge doesn't match
          // our co-kernel, something is badly wrong.
          assert(polynomialID == termID.first);
          
          const bool inserted = termIDs.insert(termID).second;
          if (!inserted)
            CFD_EXCEPTION("Duplicate term found in biclique.");

          const bool deleted = sops[polynomialID].deleteTerm(termID.second);
          if (!deleted)
            CFD_EXCEPTION("Failed to remove factorised term from SOP.");
        }
      }

      // Add new term to polynomial
      sops[polynomialID].append(get(term_cokernel(), *graph, coKernelVertex) + newCube);
    }
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
