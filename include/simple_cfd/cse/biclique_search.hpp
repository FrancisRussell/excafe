#ifndef SIMPLE_CFD_CSE_BICLIQUE_SEARCH_HPP
#define SIMPLE_CFD_CSE_BICLIQUE_SEARCH_HPP

#include <utility>
#include <cassert>
#include <algorithm>
#include <boost/foreach.hpp>
#include "biclique.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
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

  int maximumCubes;
  int maximumCubeValueSum;
  vertex_descriptor nextSplitPoint;
  bool finished;

  void calculateValues(const vertex_descriptor& oldSplitPoint)
  {
    const vertex_descriptor nullVertex = boost::graph_traits<graph_t>::null_vertex();
    int candidateCubes = 0;
    int candidateCubesValueSum = 0;
    nextSplitPoint = nullVertex;

    // nextSplitPoint is the lexicographically next cube vertex after the old split point
    // contained in the cube candidates.

    if (this->empty())
    {
      BOOST_FOREACH(const vertex_descriptor& v, vertices(this->getGraph()))
      {
        const bool afterOldSplit = (oldSplitPoint == nullVertex || v > oldSplitPoint);
        if (afterOldSplit && get(is_cube(), this->getGraph(), v))
        {
          const bool lowerThanNextSplit = (nextSplitPoint == nullVertex || v < nextSplitPoint);
          if (lowerThanNextSplit)
            nextSplitPoint = v;

          ++candidateCubes;
          candidateCubesValueSum += get(mul_count(), this->getGraph(), v);
        }
      }
    }
    else
    {
      assert(oldSplitPoint != boost::graph_traits<graph_t>::null_vertex());
      int topScore = 0;

      BOOST_FOREACH(const vertex_descriptor& coKernel, this->coKernelVertices)
      {
        int currentCandidateCubes = 0;
        int currentCandidateCubesValueSum = 0;

        BOOST_FOREACH(const edge_descriptor& e, out_edges(coKernel, this->getGraph()))
        {
          const vertex_descriptor candidateCube = target(e, this->getGraph());
          if (candidateCube > oldSplitPoint)
          {
            ++currentCandidateCubes;
            currentCandidateCubesValueSum += get(mul_count(), this->getGraph(), candidateCube);

            if (nextSplitPoint == nullVertex || candidateCube < nextSplitPoint)
              nextSplitPoint = candidateCube;
          }
        }

        const int currentScore = base_t::getValue(this->getCubeValueSum() + currentCandidateCubesValueSum, 
          this->numCubes() + currentCandidateCubes, this->getCoKernelValueSum(), this->numCoKernels());

        if (currentScore > topScore)
        {
          topScore = currentScore;
          candidateCubes = currentCandidateCubes;
          candidateCubesValueSum = currentCandidateCubesValueSum;
        }
      }
    }

    // Find next split point
    finished = (nextSplitPoint == nullVertex);

    maximumCubes = this->numCubes() + candidateCubes;
    maximumCubeValueSum = this->getCubeValueSum() + candidateCubesValueSum;
  }

  // This constructs really large sets, so keep it private.
  BicliqueSearch(graph_t& g) : base_t(g)
  {
    const vertex_descriptor nullSplitPoint = boost::graph_traits<graph_t>::null_vertex();
    calculateValues(nullSplitPoint);
  }

  BicliqueSearch(const BicliqueSearch& b, const vertex_descriptor& newVertex, const same_biclique_tag) : 
    base_t(b)
  {
    calculateValues(newVertex);
  }

  BicliqueSearch(const BicliqueSearch& b, const vertex_descriptor& newVertex, const grow_biclique_tag) :
    base_t(b, newVertex)
  {
    calculateValues(newVertex);
  }

public:
  BicliqueSearch(graph_t& graph, const vertex_descriptor& seed) : base_t(graph, seed)
  {
    calculateValues(seed);
  }

  int getMaximalValue() const
  {
    return base_t::getValue(maximumCubeValueSum, maximumCubes, this->getCoKernelValueSum(),
      this->numCoKernels());
  }

  std::pair<BicliqueSearch, BicliqueSearch> split() const
  {
    if (finished)
      CFD_EXCEPTION("Cannot split a finished BicliqueSearch.");

    return std::make_pair(BicliqueSearch(*this, nextSplitPoint, grow_biclique_tag()),
                          BicliqueSearch(*this, nextSplitPoint, same_biclique_tag()));
  }

  void print() const
  {
    std::cout << "num_cubes=" << this->numCubes() << ", num_cokernels=" << this->numCoKernels();
    std::cout << ", value=" << this->getValue() << ", maximal_value=" << getMaximalValue();
    std::cout << ", split_vertex=" << nextSplitPoint << ", candidate_cubes=" << maximumCubes;
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
    std::swap(maximumCubes, b.maximumCubes);
    std::swap(maximumCubeValueSum, b.maximumCubeValueSum);
    std::swap(nextSplitPoint, b.nextSplitPoint);
  }
};

struct BicliqueSearchComparator
{
  template<typename graph_t>
  bool operator()(const BicliqueSearch<graph_t>& a, const BicliqueSearch<graph_t>& b) const
  {
    return a.getMaximalValue() < b.getMaximalValue();
  }
};

}

}

namespace std
{

template<typename G>
void swap(cfd::cse::BicliqueSearch<G>& a, cfd::cse::BicliqueSearch<G>& b)
{
  a.swap(b);
}

}

#endif
