#ifndef SIMPLE_CFD_CSE_BICLIQUE_SEARCH_HPP
#define SIMPLE_CFD_CSE_BICLIQUE_SEARCH_HPP

#include <utility>
#include <cassert>
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
    std::set<vertex_descriptor> candidateCubes;

    // Efficiency warning: calculating the maximal value of an empty clique involves building a set
    // containing every cube vertex in the entire graph more than the old split point. This could be
    // extremely large.
    if (this->empty())
    {
      BOOST_FOREACH(const vertex_descriptor& v, vertices(this->getGraph()))
      {
        const bool isNullSplitPoint = oldSplitPoint==boost::graph_traits<graph_t>::null_vertex();
        const bool isAfterSplitPoint = isNullSplitPoint || (!isNullSplitPoint && v>oldSplitPoint);

        if (isAfterSplitPoint && get(is_cube(), this->getGraph(), v))
          candidateCubes.insert(v);
      }
    }
    else
    {
      assert(oldSplitPoint != boost::graph_traits<graph_t>::null_vertex());

      BOOST_FOREACH(const vertex_descriptor& coKernel, this->coKernelVertices)
      {
        BOOST_FOREACH(const edge_descriptor& e, out_edges(coKernel, this->getGraph()))
        {
          if (target(e, this->getGraph()) > oldSplitPoint)
          {
            candidateCubes.insert(target(e, this->getGraph()));
          }
        }
      }
    }

    // Find next split point
    if (candidateCubes.empty())
    {
      finished = true;
      nextSplitPoint = boost::graph_traits<graph_t>::null_vertex();
    }
    else
    {
      finished = false;
      nextSplitPoint = *candidateCubes.begin();
    }

    maximumCubes = this->numCubes() + candidateCubes.size();
    maximumCubeValueSum = this->getCubeValueSum() + base_t::getValue(this->getGraph(), candidateCubes.begin(),
      candidateCubes.end());
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
