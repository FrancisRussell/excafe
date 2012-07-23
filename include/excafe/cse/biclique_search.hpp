#ifndef EXCAFE_CSE_BICLIQUE_SEARCH_HPP
#define EXCAFE_CSE_BICLIQUE_SEARCH_HPP

#include <utility>
#include <cassert>
#include <algorithm>
#include <boost/foreach.hpp>
#include "biclique.hpp"
#include "vertex_info.hpp"
#include <excafe/exception.hpp>

namespace excafe
{

namespace cse
{

template<typename G>
class BicliqueSearch : public Biclique<G>
{
private:
  struct same_biclique_tag {};
  struct grow_biclique_tag {};

  typedef G graph_t;
  typedef Biclique<graph_t> base_t;
  typedef typename base_t::vertex_descriptor vertex_descriptor;
  typedef typename base_t::edge_descriptor   edge_descriptor;

  long depth;
  VertexInfo<graph_t> maxCubeInfo;
  vertex_descriptor nextSplitPoint;
  bool finished;

  typename boost::vertex_bundle_type<graph_t>::type::cube_ordering_t
  id(const vertex_descriptor& v) const
  {
    return this->getGraph()[v].cube_ordering;
  }

  void calculateValues(const vertex_descriptor& oldSplitPoint)
  {
    const vertex_descriptor nullVertex = boost::graph_traits<graph_t>::null_vertex();
    VertexInfo<graph_t> candidateCubeInfo;
    nextSplitPoint = nullVertex;

    // nextSplitPoint is the lexicographically next cube vertex after the old split point
    // contained in the cube candidates.

    if (this->empty())
    {
      BOOST_FOREACH(const vertex_descriptor& v, vertices(this->getGraph()))
      {
        const bool afterOldSplit = (oldSplitPoint == nullVertex || id(v) > id(oldSplitPoint));
        if (afterOldSplit && this->getGraph()[v].is_cube)
        {
          const bool lowerThanNextSplit = (nextSplitPoint == nullVertex || id(v) < id(nextSplitPoint));
          if (lowerThanNextSplit)
            nextSplitPoint = v;

          candidateCubeInfo.addVertex(this->getGraph(), v);
        }
      }
    }
    else
    {
      assert(oldSplitPoint != boost::graph_traits<graph_t>::null_vertex());
      int topScore = 0;

      BOOST_FOREACH(const vertex_descriptor& coKernel, *this->coKernelVertices)
      {
        VertexInfo<graph_t> currentCandidateInfo;

        BOOST_FOREACH(const edge_descriptor& e, out_edges(coKernel, this->getGraph()))
        {
          const vertex_descriptor candidateCube = target(e, this->getGraph());
          if (id(candidateCube) > id(oldSplitPoint))
          {
            currentCandidateInfo.addVertex(this->getGraph(), candidateCube);

            if (nextSplitPoint == nullVertex || id(candidateCube) < id(nextSplitPoint))
              nextSplitPoint = candidateCube;
          }
        }

        const int currentScore = base_t::getValue(this->cubeInfo + currentCandidateInfo, this->coKernelInfo);
        if (currentScore > topScore)
        {
          topScore = currentScore;
          candidateCubeInfo = currentCandidateInfo;
        }
      }
    }

    // Find next split point
    finished = (nextSplitPoint == nullVertex);
    maxCubeInfo = this->cubeInfo + candidateCubeInfo;
  }

  // This constructs really large sets, so keep it private.
  BicliqueSearch(graph_t& g) : base_t(g), depth(0)
  {
    const vertex_descriptor nullSplitPoint = boost::graph_traits<graph_t>::null_vertex();
    calculateValues(nullSplitPoint);
  }

  BicliqueSearch(const BicliqueSearch& b, const vertex_descriptor& newVertex, const same_biclique_tag) : 
    base_t(b), depth(b.depth+1)
  {
    calculateValues(newVertex);
  }

  BicliqueSearch(const BicliqueSearch& b, const vertex_descriptor& newVertex, const grow_biclique_tag) :
    base_t(b, newVertex), depth(b.depth+1)
  {
    calculateValues(newVertex);
  }

public:
  BicliqueSearch(graph_t& graph, const vertex_descriptor& seed) : base_t(graph, seed), depth(0)
  {
    calculateValues(seed);
  }

  int getMaximalValue() const
  {
    return base_t::getValue(maxCubeInfo, this->coKernelInfo);
  }

  std::pair<BicliqueSearch, BicliqueSearch> split() const
  {
    if (finished)
      CFD_EXCEPTION("Cannot split a finished BicliqueSearch.");

    // If growing this biclique by a cube doesn't reduce the number of co-kernels,
    // there is no point in considering the case where we don't add the cube.
    const BicliqueSearch grown(*this, nextSplitPoint, grow_biclique_tag());
    if(grown.numCoKernels() == this->numCoKernels() && !grown.isFinished())
    {
      return grown.split();
    }
    else
    {
      return std::make_pair(grown, BicliqueSearch(*this, nextSplitPoint, same_biclique_tag()));
    }
  }

  std::size_t getDepth() const
  {
    return depth;
  }

  void print() const
  {
    std::cout << "num_cubes=" << this->numCubes() << ", num_cokernels=" << this->numCoKernels();
    std::cout << ", value=" << this->getValue() << ", maximal_value=" << getMaximalValue();
    std::cout << ", depth=" << depth << ", split_vertex=" << nextSplitPoint;
    std::cout << ", candidate_cubes=" << maxCubeInfo.getCount();
  }

  bool isFinished() const
  {
    return finished;
  }

  Biclique<graph_t> getBiclique() const
  {
    return Biclique<graph_t>(*this);
  }

  void swap(BicliqueSearch& b)
  {
    base_t::swap(b);
    std::swap(depth, b.depth);
    std::swap(maxCubeInfo, b.maxCubeInfo);
    std::swap(nextSplitPoint, b.nextSplitPoint);
    std::swap(finished, b.finished);
  }
};

struct BicliqueSearchComparator
{
  template<typename graph_t>
  bool operator()(const BicliqueSearch<graph_t>& a, const BicliqueSearch<graph_t>& b) const
  {
    return std::make_pair(a.getMaximalValue(), -a.getDepth()) < 
           std::make_pair(b.getMaximalValue(), -b.getDepth());
  }
};

}

}

namespace std
{

template<typename G>
void swap(excafe::cse::BicliqueSearch<G>& a, excafe::cse::BicliqueSearch<G>& b)
{
  a.swap(b);
}

}

#endif
