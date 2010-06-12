#ifndef SIMPLE_CFD_CSE_BICLIQUE_SEARCH_SPACE_HPP
#define SIMPLE_CFD_CSE_BICLIQUE_SEARCH_SPACE_HPP

#include "biclique.hpp"
#include <boost/shared_ptr.hpp>
#include <set>
#include <limits>
#include <utility>
#include <vector>

namespace cfd
{

namespace cse
{

template<typename B>
class BicliqueSearchSpace
{
private:
  typedef B biclique_t;
  typedef typename biclique_t::vertex_descriptor vertex_descriptor;
  typedef typename biclique_t::edge_descriptor edge_descriptor;
  typedef typename std::vector<vertex_descriptor>::const_iterator tail_vertex_iter;

  biclique_t biclique;
  boost::shared_ptr< const std::vector<vertex_descriptor> > tailVertices;
  tail_vertex_iter tailBegin;
  tail_vertex_iter tailEnd;

  std::size_t totalCoKernelCount;
  int         totalCoKernelValue;
  int         tailVertexSum;

  void calculateCachedValues()
  {
    const typename biclique_t::graph_t& g = biclique.getGraph();
    totalCoKernelCount = 0;
    totalCoKernelValue = 0;

    std::set<vertex_descriptor> coKernelVertices;
    BOOST_FOREACH(const vertex_descriptor v, std::make_pair(tailBegin, tailBegin))
    {
      BOOST_FOREACH(const edge_descriptor e, out_edges(v, g))
      {
        coKernelVertices.insert(target(e, g));
      }
    }

    totalCoKernelCount = coKernelVertices.size();
    BOOST_FOREACH(const vertex_descriptor v, coKernelVertices)
    {
      totalCoKernelValue += get(mul_count(), g, v);
    }

    tailVertexSum = biclique_t::getValue(g, tailBegin, tailEnd);
  }

  BicliqueSearchSpace(const BicliqueSearchSpace& parent, const biclique_t& _biclique) : 
    biclique(_biclique), tailVertices(parent.tailVertices),
    tailBegin(parent.tailBegin+1), tailEnd(parent.tailEnd), 
    totalCoKernelCount(parent.totalCoKernelCount), totalCoKernelValue(parent.totalCoKernelValue),
    tailVertexSum(parent.tailVertexSum - biclique_t::getValue(biclique.getGraph(), parent.tailBegin, tailBegin))
  {
  }

public:
  template<typename InputIterator>
  BicliqueSearchSpace(const biclique_t& _biclique, const InputIterator begin, const InputIterator end) : 
    biclique(_biclique), tailVertices(new std::vector<vertex_descriptor>(begin, end)),
    tailBegin(tailVertices->begin()), tailEnd(tailVertices->end())
  {
    calculateCachedValues();
  }

  biclique_t getBiclique() const
  {
    return biclique;
  }

  int getValue() const
  {
    return biclique.getValue();
  }

  int getMaximalValue() const
  {
    const std::size_t maximalCubeCount = biclique.numCubes() + std::distance(tailBegin, tailEnd);
    const int maximalCubeValue = biclique.getCubeValueSum() + tailVertexSum;

    // In the case we have the null biclique, we use cached values for the count and sum of co-kernel values
    const std::size_t coKernelCount = biclique.empty() ? totalCoKernelCount : biclique.numCoKernels();
    const int coKernelValue = biclique.empty() ? totalCoKernelValue : biclique.getCoKernelValueSum();

    return biclique.getValue(maximalCubeValue, maximalCubeCount, coKernelValue, coKernelCount);
  }

  std::pair<BicliqueSearchSpace, BicliqueSearchSpace> split() const
  {
    if (finished())
      CFD_EXCEPTION("Search space is empty.");

    const vertex_descriptor v = *tailBegin;
    return std::make_pair(
      BicliqueSearchSpace(*this, biclique), 
      BicliqueSearchSpace(*this, biclique.addVertex(v)));
  }

  void swap(BicliqueSearchSpace& s)
  {
    std::swap(biclique, s.biclique);
    std::swap(tailVertices, s.tailVertices);
    std::swap(tailBegin, s.tailBegin);
    std::swap(tailEnd, s.tailEnd);
  }

  void print() const
  {
    biclique.print();
    std::cout << ", maximal_value=" << getMaximalValue();
    std::cout << ", tail_size: " << std::distance(tailBegin, tailEnd) << std::endl;
  }

  bool finished() const
  {
    return tailBegin == tailEnd || getMaximalValue() == 0;
  }
};

template<typename B>
struct BicliqueSearchSpaceComparator
{
  typedef B biclique_t;

  bool operator()(const BicliqueSearchSpace<biclique_t>& a, const BicliqueSearchSpace<biclique_t>& b) const
  {
    return a.getMaximalValue() < b.getMaximalValue();
  }
};

}

}

namespace std
{

template<typename B>
void swap(cfd::cse::BicliqueSearchSpace<B>& a, cfd::cse::BicliqueSearchSpace<B>& b)
{
  a.swap(b);
}

}

#endif
